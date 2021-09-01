use super::meta::{BlockMeta, Codecs, FieldMeta, FileInfo, FileMeta, FILE_INFO_SIZE};
use crate::stats::CollectStats;
use crate::{CompressTask, Compressor, SIZE_LIMIT, U32_SIZE};
use bam_tools::record::bamrawrecord::BAMRawRecord;
use bam_tools::record::fields::{
    field_type, is_data_field, var_size_field_to_index, FieldType, Fields, FIELDS_NUM,
};
use byteorder::{LittleEndian, WriteBytesExt};
use crc32fast::Hasher;
use std::borrow::{Borrow, Cow};
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{Seek, SeekFrom, Write};

pub(crate) struct BlockInfo {
    pub numitems: u32,
    pub block_size: u32,
    pub uncompr_size: usize,
    pub field: Fields,
    // Interpretation is up to the reader.
    pub max_value: Option<Vec<u8>>,
    pub min_value: Option<Vec<u8>>,
}

impl Default for BlockInfo {
    fn default() -> Self {
        Self {
            numitems: 0,
            block_size: 0,
            uncompr_size: 0,
            field: Fields::RefID,
            max_value: None,
            min_value: None,
        }
    }
}
struct BlockSink<W>
where
    W: Write + Seek,
{
    writer: W,
    compressor: Compressor,
}

impl<W> BlockSink<W>
where
    W: Write + Seek,
{
    /// Send new data for compression. Retrieve old data
    pub fn load_recv_compression(
        &mut self,
        bytes: Vec<u8>,
        block_info: BlockInfo,
    ) -> (Vec<u8>, BlockInfo) {
        // Flush already compressed data
        let compress_task = self.compressor.get_compr_block();
    }

    ///
    fn write_block(&mut self, bytes: Vec<u8>) -> std::io::Result<u64> {
        let seek_pos = self.writer.seek(SeekFrom::Current(0)).unwrap();
        self.writer.write_all(bytes).unwrap();
        Ok(seek_pos)
    }
}

struct Inner<W, T>
where
    W: Write + Seek,
{
    stats_collector: Option<Box<dyn CollectStats<T>>>,
    sink: RefCell<Sink<W>>,
    buffer: Vec<u8>,
    offset: usize,
    field: Fields,
    rec_count: u32,
}

impl<W, T> Inner<W, T>
where
    W: Write + Seek,
{
    pub fn write_data(&mut self, data: &[u8]) {
        if self.offset > 0 && self.offset + data.len() > SIZE_LIMIT {
            self.flush();
            self.offset = 0;
        }

        self.buffer[self.offset..self.offset + data.len()].clone_from_slice(data);
        self.offset += data.len();

        self.rec_count += 1;
    }

    fn flush(&mut self) {
        if self.rec_count == 0 {
            return;
        }

        self.sink.borrow_mut().write()
    }
}

trait Column {
    // Extracts and writes data from corresponding BAMRawRecord record.
    fn write_record_field(&mut self, rec: &mut BAMRawRecord) -> ();
}

/// Column containing fixed sized fields.
struct FixedColumn<W, T>(Inner<W, T>)
where
    W: Write + Seek,
    T: Clone;

impl<W, T> Column for FixedColumn<W, T>
where
    W: Write,
    T: Clone,
{
    fn write_record_field(&mut self, rec: &mut BAMRawRecord) -> () {
        let data = rec.get_bytes(&self.0.field);
        // At least one record will be written in even if it exceeds SIZE_LIMIT (for variable sized fields).

        if let Some(ref mut stats) = self.0.stats_collector {
            stats.collect_stats(rec);
        }
    }
}

/// Column containing variable sized fields.
// struct VariableColumn<W: Write> {
//     inner: Inner<W>,
//     index: FixedColumn<W>,
// }

/// The data is held in blocks.
///
/// Fixed sized fields are written as fixed size blocks into file. All blocks
/// (for fixed size fields) except last one contain equal amount of data.
///
/// Variable sized fields are written as fixed size blocks. Blocks may contain
/// different amount of data. Variable sized fields are accompanied by separate
/// index in separate block for fixed size fields. Groups records before writing
/// out to file.
pub struct Writer<W>
where
    W: Write + Seek,
{
    chunks: Vec<Vec<u8>>,
    // Current item index
    offsets: [usize; FIELDS_NUM],
    num_items: [u32; FIELDS_NUM],
    // Compression buffer to avoid allocations
    compressor: Compressor,
    file_meta: FileMeta,
    inner: W,
    // Used to order meta information, since multithreaded compressor may
    // compress latter block first, and disturb order.
    blocks_nums: Vec<usize>,
}

impl<W> Writer<W>
where
    W: Write + Seek,
{
    /// Create new writer
    pub fn new(
        mut inner: W,
        codec: Codecs,
        thread_num: usize,
        ref_seqs: Vec<(String, i32)>,
    ) -> Self {
        // Make space for the FileInfo to be written into.
        inner
            .seek(SeekFrom::Start((FILE_INFO_SIZE) as u64))
            .unwrap();
        Writer {
            chunks: vec![vec![0; SIZE_LIMIT]; FIELDS_NUM],
            offsets: [0; FIELDS_NUM],
            num_items: [0; FIELDS_NUM],
            file_meta: FileMeta::new(codec, ref_seqs),
            compressor: Compressor::new(thread_num),
            inner,
            blocks_nums: vec![0; FIELDS_NUM],
        }
    }
    /// Push BAM record into this writer
    pub fn push_record(&mut self, record: &BAMRawRecord) {
        let mut index_fields_buf: [u8; U32_SIZE] = [0; U32_SIZE];
        // Index fields are not written on their own. They hold index data for variable sized fields.
        for field in Fields::iterator().filter(|f| is_data_field(*f)) {
            let new_data = record.get_bytes(field);
            match field_type(field) {
                // Require update to index fields
                FieldType::VariableSized => {
                    // Write variable sized field
                    self.update_field_buf(field, new_data);
                    let end_pos = self.offsets[*field as usize];
                    (&mut index_fields_buf[..])
                        .write_u32::<LittleEndian>(end_pos as u32)
                        .unwrap();
                    // Write fixed size index
                    self.update_field_buf(&var_size_field_to_index(field), &index_fields_buf);
                }
                FieldType::FixedSized => {
                    self.update_field_buf(field, new_data);
                }
            }
        }
    }

    /// Used to write new data into buffers
    fn update_field_buf(&mut self, field: &Fields, new_data: &[u8]) {
        let mut offset_into_chunk = self.offsets[*field as usize];

        // At least one record will be written in even if it exceeds SIZE_LIMIT (for variable sized fields).
        if offset_into_chunk > 0 && offset_into_chunk + new_data.len() > SIZE_LIMIT {
            self.flush(field);
            offset_into_chunk = 0;
        }

        let cur_chunk = &mut self.chunks[*field as usize];
        if cur_chunk.len() < SIZE_LIMIT {
            cur_chunk.resize(std::cmp::max(new_data.len(), SIZE_LIMIT), 0);
        }
        let item_counter = &mut self.num_items[*field as usize];

        cur_chunk[offset_into_chunk..offset_into_chunk + new_data.len()].clone_from_slice(new_data);
        offset_into_chunk += new_data.len();

        self.offsets[*field as usize] = offset_into_chunk;
        *item_counter += 1;
    }

    /// This method only schedules field for compression, it doesn't immediately
    /// flush it to writer.
    fn flush(&mut self, field: &Fields) {
        // Already empty
        if self.num_items[*field as usize] == 0 {
            return;
        }

        // Flush already compressed data
        let compress_task = self.compressor.get_compr_block();

        // Skips prefilled blocks in the beginning of program execution
        if compress_task.uncompr_size != 0 {
            self.write_data_and_update_meta(&compress_task);
        }

        let mut buf = compress_task.buf;
        std::mem::swap(&mut buf, &mut self.chunks[*field as usize]);
        let uncompr_size = self.offsets[*field as usize];
        let codec = self.file_meta.get_field_codec(field);

        self.compressor.compress_block(
            self.blocks_nums[*field as usize],
            *field,
            self.num_items[*field as usize],
            uncompr_size,
            buf,
            *codec,
        );
        self.blocks_nums[*field as usize] += 1;

        self.offsets[*field as usize] = 0;
        self.num_items[*field as usize] = 0;
    }

    fn write_data_and_update_meta(&mut self, task: &CompressTask) {
        let meta = self.generate_meta(task.num_items);
        let compressed_size = task.buf.len();
        self.inner.write_all(&task.buf[..compressed_size]).unwrap();

        let block_sizes = self.file_meta.get_blocks_sizes(&task.field);
        if block_sizes.len() <= task.ordering_key {
            block_sizes.resize(task.ordering_key + 1, 0);
        }
        block_sizes[task.ordering_key] = compressed_size as u32;

        // Order as came in
        let field_meta = self.file_meta.get_blocks(&task.field);
        if field_meta.len() <= task.ordering_key {
            field_meta.resize(task.ordering_key + 1, BlockMeta::default());
        }
        field_meta[task.ordering_key] = meta;
    }

    fn generate_meta(&mut self, numitems: u32) -> BlockMeta {
        let seekpos = self.inner.seek(SeekFrom::Current(0)).unwrap();
        BlockMeta {
            seekpos,
            numitems,
            max_value: None,
            min_value: None,
        }
    }

    /// Terminates the writer. Always call after writting all the data. Returns
    /// total amount of bytes written.
    pub fn finish(&mut self) -> std::io::Result<u64> {
        // Flush leftovers
        for field in Fields::iterator() {
            self.flush(field);
        }

        for task in self.compressor.finish() {
            self.write_data_and_update_meta(&task);
        }

        let meta_start_pos = self.inner.seek(SeekFrom::Current(0))?;
        // Write meta
        let main_meta = serde_json::to_string(&self.file_meta).unwrap();
        let main_meta_bytes = main_meta.as_bytes();
        let crc32 = calc_crc_for_meta_bytes(main_meta_bytes);
        self.inner.write_all(main_meta_bytes)?;

        let total_bytes_written = self.inner.seek(SeekFrom::Current(0))?;
        // Revert back to the beginning of the file
        self.inner.seek(SeekFrom::Start(0)).unwrap();
        let file_meta = FileInfo::new([1, 0], meta_start_pos, crc32);
        let file_meta_bytes = &Into::<Vec<u8>>::into(file_meta)[..];
        self.inner.write_all(file_meta_bytes)?;
        Ok(total_bytes_written)
    }
}

impl<W> Write for Writer<W>
where
    W: Write + Seek,
{
    /// WARNING: ENSURE THAT BUF CONTAINS A ONE FULL RECORD.
    /// Write trait implementation is made to allow passing Write trait objects to sort function in BAM parallel.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        assert!(buf.len() > 0);
        // println!("Current record size is: {}", buf.len());
        let wrapper = BAMRawRecord(Cow::Borrowed(buf));
        self.push_record(&wrapper);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// TODO: Its commented out because test doesnt compile otherwise.
// impl<W> Drop for Writer<W>
// where
//     W: Write + Seek,
// {
//     fn drop(&mut self) {
//         self.finish().unwrap();
//     }
// }

pub(crate) fn calc_crc_for_meta_bytes(bytes: &[u8]) -> u32 {
    let mut hasher = Hasher::new();
    hasher.update(bytes);
    hasher.finalize()
}

#[ignore]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::parse_tmplt::*;
    use crate::reader::reader::*;
    use byteorder::ReadBytesExt;
    use std::io::Cursor;
    #[test]
    fn test_writer() {
        // let raw_records = vec![BAMRawRecord::default(); 2];
        // let mut buf: Vec<u8> = vec![0; SIZE_LIMIT];
        // let out = Cursor::new(&mut buf[..]);
        // let mut writer = Writer::new(out, Codecs::Gzip, 8);
        // for rec in raw_records.iter() {
        //     writer.push_record(rec);
        // }
        // let total_bytes_written = writer.finish().unwrap();
        // buf.resize(total_bytes_written as usize, 0);

        // let in_cursor = Box::new(Cursor::new(buf));
        // let mut parsing_template = ParsingTemplate::new();
        // parsing_template.set_all();
        // let mut reader = Reader::new(in_cursor, parsing_template).unwrap();
        // let mut records = reader.records();
        // let mut it = raw_records.iter();
        // while let Some(rec) = records.next_rec() {
        //     let rec_orig = it.next().unwrap();
        //     let orig_map_q = rec_orig.get_bytes(&Fields::Mapq)[0];
        //     let orig_pos = rec_orig
        //         .get_bytes(&Fields::Pos)
        //         .read_i32::<LittleEndian>()
        //         .unwrap();
        //     assert_eq!(rec.pos.unwrap(), orig_pos);
        //     assert_eq!(rec.mapq.unwrap(), orig_map_q);
        // }
    }
}
