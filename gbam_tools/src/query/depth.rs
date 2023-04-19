use bam_tools::record::fields::Fields;
use itertools::{Itertools, Chunk};
use rust_htslib::bam::buffer;
use std::ascii::AsciiExt;
use std::collections::BTreeMap;
use std::convert::TryInto;
use std::io::{Write, BufWriter, StdoutLock};
use std::iter::FromIterator;
use std::num;
use std::ops::{Range, RangeInclusive};
use std::sync::Arc;
use std::{cmp::Ordering, collections::HashMap, convert::TryFrom, time::Instant};
use std::fs::File;
use crate::reader::parse_tmplt::ParsingTemplate;
use crate::utils::bed;
/// This module provides function for fast querying of read depth.
use crate::meta::{BlockMeta, FileMeta};
use crate::reader::{reader::generate_block_treemap, reader::Reader, record::GbamRecord};
use std::path::PathBuf;
use rayon::prelude::*;
type Region = RangeInclusive<u32>;
use std::io::Read;
use crossbeam::channel::bounded;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::thread;
use std::thread::JoinHandle;

fn panic_err() {
    panic!("The query you entered is incorrect. The format is as following: <ref name>:<position>\ne.g. chr1:1257\n");
}
// #[derive(Clone, Default)]
// struct OperationBuffers {
//     pub increments: Vec<usize>,
//     pub decrements: Vec<usize>,
// }

fn process_range(mut gbam_reader: Reader, rec_range: RangeInclusive<usize>, mut scan_line: Vec<i32>, target_id: i32) -> Vec<i32> {
    let mut rec = GbamRecord::default();
    for idx in rec_range {
        gbam_reader.fill_record(idx, &mut rec);
        if rec.refid.unwrap() != target_id {
            continue;
        }
        let read_start: usize = rec.pos.unwrap().try_into().unwrap();
        let base_cov = rec.cigar.as_ref().unwrap().base_coverage() as usize;
        let read_end = read_start + base_cov;

        scan_line[read_start] += 1;
        scan_line[read_end] -= 1;
        // buf.increments.push(read_start);
        // buf.decrements.push(read_end);
    }
    scan_line
}

fn calc_depth(gbam_file: File, file_meta: Arc<FileMeta>, number_of_records: usize, ref_id: i32, mut coverage_arr: Vec<i32>, ref_len: usize) -> Vec<i32> {
    let lower_bound = find_leftmost_block(ref_id, file_meta.view_blocks(&Fields::RefID)).expect("RefID was not found in block meta.") as usize;
    let upper_bound = find_rightmost_block(ref_id, file_meta.view_blocks(&Fields::RefID)) as usize;
    let mut first_rec = (lower_bound as usize)*file_meta.view_blocks(&Fields::RefID)[0].numitems as usize;
    let mut last_rec = std::cmp::min(upper_bound as usize*file_meta.view_blocks(&Fields::RefID)[0].numitems as usize, number_of_records-1);

    // let mut temp_reader = Reader::new_with_meta(gbam_file.try_clone().unwrap(), ParsingTemplate::new_with(&[Fields::RefID, Fields::Pos, Fields::RawCigar]), &file_meta).unwrap();
    // let mut rec = GbamRecord::default();
    // temp_reader.fill_record(last_rec, &mut rec);

    // let read_start: usize = rec.pos.unwrap().try_into().unwrap();
    // let base_cov = rec.cigar.as_ref().unwrap().base_coverage() as usize;
    // let read_end = read_start + base_cov;

    // Loads of page faults here.
    coverage_arr.resize(ref_len, 0);

    // dbg!("Allocated {}", ref_len);

    let mut coverage = process_range(Reader::new_with_meta(gbam_file.try_clone().unwrap(), ParsingTemplate::new_with(&[Fields::RefID, Fields::Pos, Fields::RawCigar]), &file_meta).unwrap(), first_rec..=last_rec, coverage_arr, ref_id);
    let mut acc = 0;
    for slot in coverage.iter_mut() {
        acc += *slot;
        *slot = acc; 
    }
    coverage
}

pub fn main_depth(gbam_file: File, bed_file: Option<&PathBuf>, bed_cli_request: Option<String>, mapq: Option<u32>, thread_num: Option<usize>){
    let mut queries = HashMap::<String, Vec<(u32, u32)>>::new();
    if let Some(bed_path) = bed_file {
        queries = bed::parse_bed_from_file(&bed_path).expect("BED file is corrupted.");
    } 
    if let Some(query) = bed_cli_request {
        queries.extend(bed::parse_bed(&mut query.as_bytes()).unwrap().into_iter());
    }
    let qual_cutoff = mapq.unwrap_or(0);

    let mut reader = Reader::new(gbam_file.try_clone().unwrap(), ParsingTemplate::new()).unwrap();
    let file_meta = reader.file_meta.clone();
    let ref_seqs = file_meta.get_ref_seqs().clone();
    let chr_to_ref_id = get_chr_name_mapping(ref_seqs.iter().map(|(chr, _)| chr), &mut reader);
    let number_of_records = reader.amount;
    drop(reader);

    // Calculate for whole file.
    if queries.is_empty() {
        ref_seqs.iter().for_each(|(chr, len)| {queries.insert(chr.clone(), vec![(0 as u32, len-1)]);});
    }

    let longest_chr = *ref_seqs.iter().map(|(_,len)| len).max().unwrap();

    let mut buffers = vec![Vec::<i32>::new()];
    if thread_num.is_some(){
        buffers = vec![Vec::<i32>::new();std::cmp::min(thread_num.unwrap(), 8)];
    }


    let mut circular_buf_channels: Vec::<Option<JoinHandle<(String, Vec<i32>)>>> = Vec::new();
    (0..buffers.len()).for_each(|_|circular_buf_channels.push(None));

    let mut idx = 0;
    let mut coverage_arr: Vec<i64> = Vec::new(); 
    coverage_arr.reserve(longest_chr as usize);

    let mut printer = ConsolePrinter::new();
    let mut iter = ref_seqs.iter();
    let mut accum = 0;     
    
    loop {
        // dbg!(buffers.len()); 
        if idx == circular_buf_channels.len() {
            idx = 0;
        }
        if circular_buf_channels[idx].is_some() {
            let (thread_chr, mut coverage_arr) = circular_buf_channels[idx].take().unwrap().join().unwrap();

            if let Some(bed_regions) = queries.get(&thread_chr) {
                // coverage_arr.resize(*ref_len as usize, 0);
                // let ref_id = chr_to_ref_id.get(chr).unwrap().unwrap();
                // buffers = calc_depth(gbam_file.try_clone().unwrap(), file_meta.clone(), number_of_records, ref_id, &mut coverage_arr, buffers);
    

                let now = Instant::now();

    
                printer.set_chr(thread_chr.clone());
               
                for bed_region in bed_regions {
                    for coord in bed_region.0..=bed_region.1 {
                        if coverage_arr[coord as usize] > 0 {
                            printer.write_efficient(coord as u64, coverage_arr[coord as usize] as i64);
                        }
                    }
                }
                accum += now.elapsed().as_millis();
    
                coverage_arr.clear();
            }
            
            buffers.push(coverage_arr);
        }

        let next_chr = iter.next(); 
        
        if let Some((chr, ref_len)) = next_chr {
            let ref_id = chr_to_ref_id.get(chr).unwrap().unwrap();
            let buf = buffers.pop().unwrap();
            let meta = file_meta.clone();
            let file = gbam_file.try_clone().unwrap();
            let t_chr = chr.clone();
            let t_ref_len = *ref_len as usize;
            let handle = thread::spawn(move || {
                (t_chr, calc_depth(file, meta, number_of_records, ref_id, buf, t_ref_len))
            });
    
            circular_buf_channels[idx] = Some(handle);
        }

        if buffers.len() == circular_buf_channels.len(){
            break;
        }

        idx += 1;
    }

    circular_buf_channels.clear();

    dbg!(accum);
    // Shouldn't allocate more.
    assert!(coverage_arr.capacity() == longest_chr as usize);
}

// TODO: Merge bed regions with tolerance to form super regions. Then do
// sweepline for each of the super regions. Print the depth results for each of
// the regions by taking chunks from the super regions or doing in process while
// calculating to avoid double iterating the super region sweep line array.

// fn parse_bytes(bytes: &Vec<u8>) -> i32 {
//     (&bytes[..]).read_i32::<LittleEndian>().unwrap()
// }

fn get_chr_name_mapping<'a, I>(ref_ids: I, reader: &mut Reader) -> HashMap<String, Option<i32>>
where
    I: Iterator<Item = &'a String>,
{
    let mut name_to_ref_id = HashMap::<String, i32>::new();
    reader
        .file_meta
        .get_ref_seqs()
        .into_iter()
        .enumerate()
        .for_each(|(i, (name, _))| {
            name_to_ref_id.insert(name.clone(), i as i32);
        });
    // https://github.com/biod/sambamba/blob/3eff9a2d8bb3097b92c72752be3c6b42dd1c59b7/BioD/bio/std/hts/bam/read.d#L902
    ref_ids
        .map(|id| (id.to_owned(), name_to_ref_id.get(id).cloned()))
        .collect::<HashMap<String, Option<i32>>>()
}

// Union of bed regions with tolerance.
// struct SuperRegion {
//     region: Region,
//     bed_regions: Vec<Region>,
// }

// Creates super regions from multiple bed regions, if they are close enough
// (within tolerance). Later the array of size of this super region will be used
// to calculate depth for each one of the nested bed regions.
// fn merge_regions(
//     regions: &Vec<(String, u32, u32)>,
//     tolerance: u32,
// ) -> HashMap<String, Vec<SuperRegion>> {
//     let ref_id_groups = regions
//         .iter()
//         .map(|a| (a.0.clone(), (a.1..=a.2)))
//         .into_iter()
//         .into_group_map();
//     let mut ret = HashMap::<String, Vec<SuperRegion>>::new();
//     for (ref_id, mut bed_regions) in ref_id_groups.into_iter() {
//         bed_regions.sort_by(|a, b| a.start().cmp(b.start()));
//         let mut consumed_regions = vec![bed_regions.first().unwrap().clone()];
//         let mut super_start = *bed_regions.first().unwrap().start();
//         let mut super_end = *bed_regions.first().unwrap().end();
//         for range in bed_regions.into_iter().skip(1) {
//             if *range.start() > super_end + tolerance {
//                 ret.entry(ref_id.clone()).or_default().push(SuperRegion {
//                     region: super_start..=super_end,
//                     bed_regions: consumed_regions,
//                 });
//                 consumed_regions = Vec::<Region>::new();
//                 super_start = *range.start();
//             }
//             super_end = std::cmp::max(super_end, *range.end());
//             consumed_regions.push(range);
//         }
//         ret.entry(ref_id.clone()).or_default().push(SuperRegion {
//             region: super_start..=super_end,
//             bed_regions: consumed_regions,
//         });
//     }
//     ret
// }

// fn get_refid_bounds(
//     mut ref_ids: Vec<i32>,
//     reader: &mut Reader,
// ) -> HashMap<i32, Option<Range<usize>>> {
//     ref_ids.sort();
//     let mut refid_bounds = HashMap::<i32, Option<Range<usize>>>::new();
//     let meta_blocks = reader.file_meta.view_blocks(&Fields::RefID);
//     // If stats were collected for this file it's possible to narrow binary search.
//     if meta_blocks.first().unwrap().stats.is_some() {
//         dbg!(meta_blocks.len());
//         let mut meta_iter = meta_blocks.iter();
//         let mut cur_offset: usize = 0;
//         for ref_id in ref_ids.into_iter() {
//             // Iterate thorugh blocks, until one with max value bigger than cur_chr is found.
//             while let Some(meta_block) = meta_iter.next() {
//                 if &meta_block.stats.as_ref().unwrap().max_value >= &ref_id {
//                     if &meta_block.stats.as_ref().unwrap().min_value > &ref_id {
//                         refid_bounds.insert(ref_id, None);
//                     } else {
//                         // This block may contain the REFID. Later on we will search this block to determine this ultimately.
//                         refid_bounds.insert(
//                             ref_id,
//                             Some(cur_offset..(cur_offset + meta_block.numitems as usize)),
//                         );
//                     }
//                     break;
//                 }
//                 cur_offset += meta_block.numitems as usize;
//             }
//         }
//     }
//     refid_bounds
// }

pub trait DepthWrite {
    fn write_depth(&self, chr: &str, coord: u64, depth: i64);
}

// fn process_depth_query<W: DepthWrite>(
//     reader: &mut Reader,
//     output: &mut W,
//     buf: &mut Vec<i32>,
//     super_regions: Vec<SuperRegion>,
//     first_record: usize,
//     ref_id: i32,
// ) {
//     let mut counter = 0;
//     let mut gbam_buffer_record = GbamRecord::default();
//     gbam_buffer_record.cigar = Some(super::cigar::Cigar(Vec::with_capacity(100000)));
//     for SuperRegion {
//         region: super_region,
//         bed_regions,
//     } in super_regions.into_iter()
//     {
//         counter += 1;
//         // dbg!(counter);
//         buf.clear();
//         buf.resize((super_region.end() - super_region.start() + 1) as usize, 0);
//         calc_depth(
//             reader,
//             first_record,
//             ref_id,
//             &super_region,
//             buf,
//             &mut gbam_buffer_record,
//         );
//         let offset = *super_region.start();
//         for bed_region in bed_regions.into_iter() {
//             for pos in bed_region {
//                 debug_assert!(pos >= offset);
//                 output.write_depth(&(ref_id, pos, buf[(pos - offset) as usize] as u32));
//             }
//         }
//     }
// }

struct ConsolePrinter<'a> {
    buffer: [u8; 400],
    cur_chr: String,
    stdout: BufWriter<StdoutLock<'a>>,
}
impl<'a> ConsolePrinter<'a> {
    pub fn new() -> Self {
        let stdout = std::io::stdout();
        let stdout = stdout.lock();
        let stdout = BufWriter::with_capacity(32 * 1024, stdout);
        Self {  
            buffer: [0;400],
            cur_chr: "ERROR IN PROGRAM".to_owned(),
            stdout
        }
    }

    /// Saves chr name into internal buffer to achieve speedup on conversion.
    pub fn set_chr(&mut self, chr: String) {
        self.cur_chr = chr.chars().rev().collect();
    }

    pub fn write_efficient(&mut self, mut coord: u64, mut depth: i64){
        let mut ptr = self.buffer.len()-1;

        self.buffer[ptr] = '\n' as u8;
        ptr -= 1;
        while depth > 0 {
            self.buffer[ptr] = '0' as u8+(depth%10) as u8;
            depth/=10;
            ptr -= 1;
        }
        self.buffer[ptr] = '\t' as u8;
        ptr -= 1;
        while coord > 0 {
            self.buffer[ptr] = '0' as u8+(coord%10) as u8;
            coord/=10;
            ptr -= 1;
        }
        self.buffer[ptr] = '\t' as u8;
        ptr -= 1;

        for &ch in self.cur_chr.as_bytes() {
            self.buffer[ptr] = ch;
            ptr-=1;
        }

        self.stdout.write_all(&self.buffer[(ptr+1)..]).unwrap();
    }
}
impl<'a> DepthWrite for ConsolePrinter<'a> {
    fn write_depth(&self, chr: &str, coord: u64, depth: i64) {
        println!(
            "{:?}\t{}\t{}",
            chr,
            coord,
            depth
        );
    }
}

// pub fn get_regions_depths(reader: &mut Reader, regions: &Vec<(String, u32, u32)>) {
//     let ref_id_to_chr = reader
//     .file_meta
//     .get_ref_seqs()
//     .iter().enumerate()
//     .map(|(k, (refid, _))| (k as i32, refid.clone()))
//     .collect();
//     let mut printer = ConsolePrinter::new(ref_id_to_chr);
//     get_regions_depths_with_printer(reader, regions, &mut printer);
// }
// chrM    15268   438
// chrM    15269   439
// chrM    15270   425
// chrM    15271   420
// chrM    15272   418
// chrM    15273   407
// chrM    15274   268
// chrM    15275   278
// chrM    15276   281
// Approach as in https://github.com/brentp/mosdepth
/// Get depth at position. The file should be sorted!
// pub fn get_regions_depths_with_printer<W: DepthWrite>(
//     reader: &mut Reader,
//     regions: &Vec<(String, u32, u32)>,
//     output: &mut W,
// ) {
//     if !reader
//         .parsing_template
//         .check_if_active(&[Fields::RefID, Fields::Pos, Fields::RawCigar])
//     {
//         panic!("The reader should have parsing template which includes REFID, POS and CIGAR.");
//     }

//     let super_regions = merge_regions(regions, 300_000);
//     let ref_id_to_chr = get_chr_name_mapping(super_regions.keys(), reader);
//     let depth_queries = super_regions
//         .into_iter()
//         .map(|(k, v)| (ref_id_to_chr.get(&k).unwrap().unwrap(), v))
//         .collect::<BTreeMap<i32, Vec<SuperRegion>>>();

//     let ref_ids_bounds =
//         get_refid_bounds(ref_id_to_chr.values().copied().flatten().collect(), reader);

//     let mut buf = Vec::<i32>::new();
//     for (ref_id, super_regions) in depth_queries.into_iter() {
//         if let Some(ref_id_pos_hint) = ref_ids_bounds.get(&ref_id).unwrap() {
//             // Find number of record when records with this ref_id begin.
//             // We don't need to fetch anything except refid to find first record with requested refid, so disable other fields fetching.
//             reader.fetch_only(&[Fields::RefID]);
//             let first_record = find_refid(reader, ref_id, ref_id_pos_hint);
//             reader.restore_template();
//             if let Some(first_pos) = first_record {
//                 process_depth_query(reader, output, &mut buf, super_regions, first_pos, ref_id);
//             }
//         }
//     }
// }

fn find_leftmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> Option<i64> {
    let mut left: i64 = -1;
    let mut right: i64 = block_metas.len() as i64 +1;
    while (right-left) > 1 {
        let mid = (left + right) / 2;
        let max_val = &block_metas[mid as usize].stats.as_ref().unwrap().max_value;
        match max_val.cmp(&id) {
            Ordering::Equal | Ordering::Greater => right = mid,
            Ordering::Less => left = mid,
        }
    }
    if block_metas[right as usize].stats.as_ref().unwrap().min_value > id {
        return None;
    }
    Some(right)
}

fn find_rightmost_block(id: i32, block_metas: &Vec<BlockMeta>) -> i64 {
    let mut left: i64 = -1;
    let mut right: i64 = block_metas.len() as i64;
    while (right-left) > 1 {
        let mid = (left + right) / 2;
        let min_val = &block_metas[mid as usize].stats.as_ref().unwrap().min_value;
        match min_val.cmp(&id) {
            Ordering::Equal | Ordering::Less => left = mid,
            Ordering::Greater => right = mid,
        }
    }
    right
}



// For each guess there may be I/O operation with decompression, so this method is not fast.
// fn find_refid(reader: &mut Reader, chr_num: i32, range: &Range<usize>) -> Option<usize> {
//     let pred = |num: usize, buf: &mut GbamRecord| {
//         reader.fill_record(num, buf);
//         buf.refid.unwrap().cmp(&chr_num)
//     };
//     // dbg!("Looking for", chr_num);
//     binary_search(range.start, range.end, pred)
// }

// Searches for the first record which satisfies predicate.
// fn binary_search<F: FnMut(usize, &mut GbamRecord) -> Ordering>(
//     mut left: usize,
//     mut right: usize,
//     mut cmp: F,
// ) -> Option<usize> {
//     let mut buf = GbamRecord::default();
//     let end = right;

//     while left < right {
//         let mid = (left + right) / 2;
//         // println!("Mid - {:?}", mid);
//         // dbg!(mid, &buf.refid);
//         match cmp(mid, &mut buf) {
//             Ordering::Less => left = mid + 1,
//             Ordering::Equal | Ordering::Greater => right = mid,
//         }
//     }
//     // dbg!("Finished search");
//     if left == end || cmp(left, &mut buf) != Ordering::Equal {
//         return None;
//     }
//     Some(left)
// }

// pub fn calc_depth(
//     reader: &mut Reader,
//     mut record_num: usize,
//     refid: i32,
//     super_region: &RangeInclusive<u32>,
//     buffer: &mut Vec<i32>,
//     buf: &mut GbamRecord,
// ) {
//     let sweeping_line = buffer;
//     let mut cur_depth = 0;
//     // dbg!("What?");
//     loop {
//         if record_num >= reader.amount {
//             break;
//         }

//         reader.fill_record(record_num, buf);
//         record_num += 1;

//         let read_start = buf.pos.unwrap().try_into().unwrap();
//         let base_cov = buf.cigar.as_ref().unwrap().base_coverage();
//         let read_end = read_start + base_cov;

//         // println!("That is test! {} end {}", read_start, read_end);

//         if read_start >= *super_region.end() || buf.refid.unwrap() > refid {
//             break;
//         }

//         if super_region.contains(&read_start) {
//             sweeping_line[(read_start - super_region.start()) as usize] += 1;
//         }
//         // Adding 1 since read_end is also covered.
//         if super_region.contains(&(read_end + 1)) {
//             // Record starts before the region and end on it.
//             if !super_region.contains(&read_start) {
//                 cur_depth += 1;
//             }
//             sweeping_line[(read_end + 1 - super_region.start()) as usize] -= 1;
//         }
//     }
//     // dbg!("What?");
//     let depth_of_region = sweeping_line;

//     for sweep in depth_of_region.iter_mut() {
//         debug_assert!(!(*sweep < 0 && -*sweep > cur_depth));
//         cur_depth += *sweep;
//         *sweep = cur_depth;
//     }
// }

// println!(
//     "REF_ID: {}\nPOS: {}\nbase_cov: {}\ncigar_len: {}\nsearch pos: {}\nEND: {}\nDEPTH: {}\n",
//     buf.refid.unwrap(),
//     read_pos,
//     base_cov,
//     buf.cigar.as_ref().unwrap().0.len(),
//     pos,
//     read_pos + i32::try_from(base_cov).unwrap(),
//     read_pos > pos
// );
