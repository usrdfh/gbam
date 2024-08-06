use std::collections::BTreeMap;
use std::sync::Arc;
use std::{borrow::Borrow, fs::File};

use bam_tools::record::fields::{
    field_type, var_size_field_to_index, FieldType, Fields, FIELDS_NUM,
};
use memmap2::MmapOptions;
use memmap2::Mmap;
use rust_htslib::htslib::bam1_t;

use crate::meta::{FileInfo, FileMeta, FILE_INFO_SIZE, BlockMeta};
use crate::writer::calc_crc_for_meta_bytes;

use super::{
    column::{Column, FixedColumn, Inner, VariableColumn},
    parse_tmplt::ParsingTemplate,
    record::GbamRecord,
    record::gen_hts_rec,
    records::Records,
};

use std::convert::TryFrom;

pub struct Reader {
    // Instead of hashmap. Empty columns will contain None.
    pub columns: Vec<Option<Box<dyn Column + Send>>>,
    pub parsing_template: ParsingTemplate,
    original_template: ParsingTemplate,
    pub amount: usize,
    pub next_not_read: usize,
    pub file_meta: Arc<FileMeta>,
    // Kept so File won't drop while used by mmap.
    _inner: Box<File>,
    index_mapping: Option<Arc<Vec<u32>>>,
    pub mmap: Arc<Mmap>
}

impl Reader {
    pub fn new(inner: File, parsing_template: ParsingTemplate) -> std::io::Result<Self> {
        let inner = inner;
        let mmap = unsafe { Mmap::map(inner.borrow())? };
        let file_meta = verify_and_parse_meta(&mmap)?;
        Self::new_with_meta(inner, parsing_template, &Arc::new(file_meta), None)
    }

    pub fn new_with_index(inner: File, parsing_template: ParsingTemplate, index_mapping: Option<Arc<Vec<u32>>>) -> std::io::Result<Self> {
        let inner = inner;
        let mmap = unsafe { Mmap::map(inner.borrow())? };
        let file_meta = verify_and_parse_meta(&mmap)?;
        Self::new_with_meta(inner, parsing_template, &Arc::new(file_meta), index_mapping)
    }

    pub fn new_with_meta(_inner: File, parsing_template: ParsingTemplate, file_meta: &Arc<FileMeta>, index_mapping: Option<Arc<Vec<u32>>>) -> std::io::Result<Self> {
        let _copy = _inner.try_clone()?;
        let _inner: Box<File> = Box::new(_inner);
        
        let mmap = Arc::new(unsafe { MmapOptions::new().map(&_copy)? });
        // mmap.advise(memmap2::Advice::WillNeed)?;
        // Consumes up to 16 percent of runtime on big files (20GB).
        // verify(&mmap)?;
        let amount = usize::try_from(file_meta
            .view_blocks(&Fields::RefID)
            .iter()
            .fold(0, |acc: u64, x| acc + u64::from(x.numitems))).unwrap();
        let meta = file_meta.clone();

        
        Ok(Self {
            columns: init_columns(&mmap, &parsing_template, &meta),
            original_template: parsing_template.clone(),
            parsing_template,
            file_meta: meta,
            amount,
            _inner,
            mmap,
            index_mapping: index_mapping.clone(),
            next_not_read: 0,
        })
    }

    #[inline(always)]
    pub fn fill_record(&mut self, mut rec_num: usize, rec: &mut GbamRecord) {
        // That's for index sorting
        if let Some(index_map) = &self.index_mapping {
            rec_num = index_map[rec_num] as usize;
        }
        assert!(rec_num < self.amount);
        for &field in self.parsing_template.get_active_data_fields_iter() {
            self.columns[field as usize]
                .as_mut()
                .unwrap()
                .fill_record_field(rec_num, rec);
        }
    }

    #[inline(always)]
    pub fn fill_bam1_t_record(&mut self, mut rec_num: usize, rec: &mut bam1_t) {
        // That's for index sorting
        if let Some(index_map) = &self.index_mapping {
            rec_num = index_map[rec_num] as usize;
        }
        assert!(rec_num < self.amount);

        let mut set_fields = vec![false; FIELDS_NUM];
        for r in self.parsing_template.get_active_data_fields_iter() {
            set_fields[*r as usize] = true;
        }
        let mut temp : Vec<Option::<&[u8]>> = vec![None; FIELDS_NUM];
        for (i, m) in self.columns.iter_mut().enumerate() {
            if set_fields[i] {
                temp[i] = Some(m.as_mut().unwrap().get_field_bytes(rec_num));
            }
        }

        gen_hts_rec(temp, rec)
    }

    pub fn get_column(&mut self, field: &Fields) -> &mut Box<dyn Column + Send> {
        self.columns[*field as usize]
            .as_mut()
            .unwrap()
    }

    // Temporarily disable fetching for fields which are not needed
    pub fn fetch_only(&mut self, fields: &[Fields]) {
        self.parsing_template.clear();
        for field in fields {
            self.parsing_template.set(field, true);
        }
    }

    // Restores original template if some fields fetching was paused.
    pub fn restore_template(&mut self) {
        self.parsing_template = self.original_template.clone();
    }

    /// Get iterator over all GBAM records (according to parsing template).
    pub fn records(&mut self) -> Records {
        Records::new(self)
    }
}

fn init_columns(
    mmap: &Arc<Mmap>,
    parse_template: &ParsingTemplate,
    meta: &Arc<FileMeta>,
) -> Vec<Option<Box<dyn Column + Send>>> {
    let mut res = Vec::new();
    (0..FIELDS_NUM).for_each(|_| res.push(None));
    for &field in parse_template.get_active_fields_iter() {
        res[field as usize] = Some(init_col(field, mmap, meta));
    }
    res
}

fn init_col(field: Fields, mmap: &Arc<Mmap>, meta: &Arc<FileMeta>) -> Box<dyn Column + Send> {
    let inner = Inner::new(meta.clone(), field, mmap.clone());
    match field_type(&field) {
        FieldType::FixedSized => Box::new(FixedColumn::new(inner, meta.get_field_size(&field).unwrap() as usize)),
        FieldType::VariableSized => {
            let idx_field = var_size_field_to_index(&field);
            let idx_inner = Inner::new(meta.clone(), idx_field, mmap.clone());
            let idx_col = FixedColumn::new(idx_inner, meta.get_field_size(&idx_field).unwrap() as usize);
            Box::new(VariableColumn::new(inner, idx_col))
        }
    }
}

fn parse_file_info(mmap: &Mmap) -> FileInfo {
    let file_info_bytes = &mmap[0..FILE_INFO_SIZE];
    let end_of_json = file_info_bytes.iter().position(|&r| r == 0).unwrap();
    let file_info_str = String::from_utf8(file_info_bytes[..end_of_json].to_owned()).unwrap();
    serde_json::from_str(&file_info_str).expect("File meta json string was damaged.")
}

#[allow(dead_code)]
fn verify(mmap: &Mmap) -> std::io::Result<()>{
    let file_info = parse_file_info(mmap);
    // Read file meta
    let buf = &mmap[file_info.seekpos as usize..];
    if calc_crc_for_meta_bytes(buf) != file_info.crc32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Metadata JSON was damaged.",
        ));
    }
    Ok(())
}
fn verify_and_parse_meta(mmap: &Mmap) -> std::io::Result<FileMeta> {
    let file_info = parse_file_info(mmap);
    // Read file meta
    let buf = &mmap[file_info.seekpos as usize..];
    if calc_crc_for_meta_bytes(buf) != file_info.crc32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Metadata JSON was damaged.",
        ));
    }
    let file_meta_json_str = String::from_utf8(buf.to_owned()).unwrap();
    Ok(serde_json::from_str(&file_meta_json_str).expect("File meta json string was damaged."))
}

// The tree map will be used to quickly determine which block record belong to.
pub(crate) fn generate_block_treemap(meta: &FileMeta, field: &Fields) -> BTreeMap<usize, usize> {
    meta.view_blocks(field)
        .iter()
        .enumerate()
        // Prefix sum.
        .scan(0, |acc: &mut u64, (block_index, x): (usize, &BlockMeta)| {
            let current_chunk = Some((usize::try_from(*acc).unwrap() , block_index));
            *acc += x.numitems as u64;
            current_chunk
        })
        .collect()
}
