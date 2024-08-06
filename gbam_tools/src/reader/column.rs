use std::{collections::BTreeMap, io::Result, sync::Arc};

use super::reader::generate_block_treemap;
use super::record::GbamRecord;
use crate::SIZE_LIMIT;
use lzzzz::{lz4};
use bam_tools::record::fields::Fields;
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::write::GzDecoder;
use memmap2::Mmap;
use std::convert::TryFrom;

use crate::{meta::FileMeta, Codecs};

// Contains fields needed both for fixed sized fields and variable sized fields.
pub struct Inner {
    /// Arc is needed since this struct should work with PyO3 which sends struct between threads (Send trait is required).
    meta: Arc<FileMeta>,
    range_begin: usize,
    range_end: usize,
    field: Fields,
    buffer: Vec<u8>,
    reader: Arc<Mmap>,
}

impl Inner {
    pub(crate) fn new(meta: Arc<FileMeta>, field: Fields, reader: Arc<Mmap>) -> Self {
        Inner {
            meta,
            range_begin: 0,
            range_end: 0,
            field,
            buffer: Vec::<u8>::with_capacity(SIZE_LIMIT * 2),
            reader,
        }
    }
}

/// Defines how columns will operate. It is needed since variable sized fields
/// columns also require parsing of additional fixed sized fields columns.
pub trait Column {
    // Fills GbamRecord field with data from corresponding BAM record.
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) ;
    // For more direct coordination of calls.
    fn get_field_bytes(&mut self, item_num: usize) -> &[u8];
}

/// GBAM file column. Responsible for fetching data.
pub struct FixedColumn(Inner, usize);

impl Column for FixedColumn {
    /// Fetches data into provider record buffer. If item is located outside of
    /// currently loaded data block, the new block will be loaded and
    /// decompressed.
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) {
        rec.parse_from_bytes(&self.0.field.clone(), self.get_item(item_num));
    }

    fn get_field_bytes(&mut self, item_num: usize) -> &[u8] {
        self.get_item(item_num)
    }
}

impl FixedColumn {
    pub fn new(inner: Inner, field_size: usize) -> Self {
        Self(inner, field_size)
    }
    fn get_item(&mut self, item_num: usize) -> &[u8] {
        if let Some(block_num) = self.find_block(item_num) {
            Self::update_buffer(&mut self.0, block_num);
        }
        let rec_num_in_block = item_num - self.0.range_begin;
        let item_size = self.1;
        let offset = rec_num_in_block * item_size;
        &self.0.buffer[offset..offset + item_size]
    }
    // Finds blocks where record is located. None is returned if block is already loaded.
    fn find_block(&self, item_num: usize) -> Option<usize> {
        if item_num >= self.0.range_begin && item_num < self.0.range_end {
            return None;
        }
        // All blocks sizes are equal except maybe the last one since it's a fixed sized column and block size limit is constant.
        let block_len = self.0.meta.view_blocks(&self.0.field)[0].numitems;
        Some(item_num / block_len as usize)
    }

    fn update_buffer(inner: &mut Inner, block_num: usize) {
        fetch_block(inner, block_num).unwrap();
        let block_len = inner.meta.view_blocks(&inner.field)[0].numitems as usize;
        let cur_block_len = inner.meta.view_blocks(&inner.field)[block_num].numitems as usize;
        inner.range_begin = block_num * block_len;
        inner.range_end = inner.range_begin + cur_block_len;
    }
}

/// Column managing access to variable sized data. Utilizes another column (for fixed sized fields) to index data.
pub struct VariableColumn {
    inner: Inner,
    index: FixedColumn,
    // Used to quickly determine what block record belongs to.
    blocks: BTreeMap<usize, usize>,
}

impl Column for VariableColumn {
    fn fill_record_field(&mut self, item_num: usize, rec: &mut GbamRecord) {
        rec.parse_from_bytes(&self.inner.field.clone(), self.get_item(item_num));
    }

    fn get_field_bytes(&mut self, item_num: usize) -> &[u8] {
        self.get_item(item_num)
    }
}

impl VariableColumn {
    pub fn new(inner: Inner, index: FixedColumn) -> Self {
        Self {
            blocks: generate_block_treemap(&inner.meta, &inner.field),
            inner,
            index,
        }
    }

    fn get_item(&mut self, item_num: usize) -> &[u8] {
        if let Some((range_begin, block_num)) = self.find_block(item_num) {
            Self::update_buffer(&mut self.inner, block_num, range_begin);
        }
        let rec_num_in_block = item_num - self.inner.range_begin;
        let mut read_offset =
            |n| self.index.get_item(n).read_u32::<LittleEndian>().unwrap() as usize;
        let start = match rec_num_in_block {
            0 => 0,
            _ => read_offset(item_num - 1),
        };
        let end = read_offset(item_num);
        &self.inner.buffer[start..end]
    }

    // Finds blocks where record is located. None is returned if block is already loaded.
    fn find_block(&self, item_num: usize) -> Option<(usize, usize)> {
        if item_num >= self.inner.range_begin && item_num < self.inner.range_end {
            return None;
        }
        // To determine what block record N is in.
        Some(
            self.blocks
                // Inclusive range.
                .range(..=item_num)
                .next_back()
                .map_or((0, 0), |(&range_begin, &block_num)| {
                    (range_begin, block_num)
                }),
        )
    }

    fn update_buffer(inner: &mut Inner, block_num: usize, range_begin: usize) {
        fetch_block(inner, block_num).unwrap();
        let block_len = inner.meta.view_blocks(&inner.field)[block_num].numitems as usize;
        inner.range_begin = range_begin;
        inner.range_end = inner.range_begin + block_len;
    }
}

/// Fetch and decompress a data block.
fn fetch_block(inner_column: &mut Inner, block_num: usize) -> Result<()> {
    // println!("Fetching for {}", inner_column.field);
    let field = &inner_column.field;
    let block_meta = inner_column.meta.view_blocks(field).get(block_num).unwrap();
    let reader = &inner_column.reader;
    let block_size = block_meta.block_size;
    let uncompressed_size = block_meta.uncompressed_size;

    let data =
        &reader[usize::try_from(block_meta.seekpos).unwrap()..usize::try_from(block_meta.seekpos + block_size as u64).unwrap()];
    // inner_column.buffer.clear();
    // dbg!(uncompressed_size);
    inner_column.buffer.resize(uncompressed_size as usize, 0);
    let codec = inner_column.meta.get_field_codec(field);

    if uncompressed_size > 0 {
        decompress_block(data, &mut inner_column.buffer, codec).expect("Decompression failed.");
    }
    
    Ok(())
}


pub fn decompress_block(source: &[u8], dest: &mut Vec<u8>, codec: &Codecs) -> std::io::Result<()> {
    use std::io::Write;
    match codec {
        Codecs::Gzip => {
            let mut decoder = GzDecoder::new(dest);
            decoder.write_all(source).unwrap();
            decoder.try_finish().unwrap();
        }
        Codecs::Lz4 => {
            lz4::decompress(source, dest).unwrap();
        }
        Codecs::NoCompression => {
            dest.clear();
            dest.extend_from_slice(source);
        }
    };
    Ok(())
}
