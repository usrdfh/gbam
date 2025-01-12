// use gbam_tools::bam_to_gbam;
use bam_tools::{record::fields::Fields, MEGA_BYTE_SIZE};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use gbam_tools::{
    bam::bam_to_gbam::bam_sort_to_gbam,
    bam::gbam_to_bam::gbam_to_bam,
    query::depth::main_depth,
    reader::{parse_tmplt::ParsingTemplate, reader::Reader, record::GbamRecord},
    {bam_to_gbam, Codecs},
    query::flagstat::collect_stats,
};
use itertools::zip_eq;
use std::fs::OpenOptions;

use std::io::{BufRead, BufReader, Seek, SeekFrom};


use std::{path::PathBuf, convert::TryInto, io::{Read}, io::{BufWriter, Write}};
use std::time::Instant;
use std::fs::File;
use structopt::StructOpt;
use std::env;

use gbam_tools::query::cigar::base_coverage;

use rayon::prelude::*;

#[derive(StructOpt)]
struct Cli {
    /// Sort BAM file before converting it to GBAM.
    #[structopt(short, long)]
    sort: bool,
    /// Specify which kind of temporary medium to use while sorting: ram, lz4_ram, file, lz4_file
    #[structopt(long)]
    sort_temp_mode: Option<String>,
    /// Determines whether conversion is requested
    #[structopt(short, long)]
    convert_to_gbam: bool,
    /// Convert to bam
    #[structopt(long)]
    convert_to_bam: bool,
    /// Perform the test
    #[structopt(short, long)]
    test: bool,
    /// Fetch cigar in parallel for testing purposes.
    #[structopt(short, long)]
    parallel_cigar_fetch: bool,
    /// Get depth at position.
    #[structopt(short, long)]
    depth: bool,
    /// Collect statistic from flag field from all records in the file.
    #[structopt(short, long)]
    flagstat: bool,
    /// The path to the BAM file to read
    #[structopt(parse(from_os_str))]
    in_path: PathBuf,
    /// The path to write output GBAM file
    #[structopt(short, parse(from_os_str))]
    out_path: Option<PathBuf>,
    /// Depth query. Example: chr1:54-54, or chrX:1258-9999
    #[structopt(short, long)]
    query: Option<String>,
    /// Depth query. Example: chr1:54, or chrX:1258
    #[structopt(short, parse(from_os_str))]
    bed_file: Option<PathBuf>,
    /// Depth query. Filter reads with map quality lower than.
    #[structopt(long)]
    mapq: Option<u32>,
    /// Depth query. Number of threads to use. WARNING: each thread will attempt to allocate up to 1GB.
    #[structopt(long)]
    thread_num: Option<usize>,
    /// Sort temp directory.
    #[structopt(long, parse(from_os_str))]
    temp_dir: Option<PathBuf>,
    /// View header
    #[structopt(short, long)]
    header: bool,
    /// View file in binary format. Can be piped to samtools view. `gbam_binary -v test_data/1gb.gbam | samtools view`
    #[structopt(short, long)]
    view: bool,
    /// View file in binary format for piping to samtools markdup. `gbam_binary -v little.gbam > /tmp/testpipe.bam & samtools markdup -u /tmp/testpipe.bam /tmp/testoutpipe.bam`. It disables reading of two heavy fields to potentially speedup the process.
    #[structopt(long)]
    markdup_view: bool,
    /// Patch the gbam file with dups detected by other software. A multiline file is expected with 0 and 1, where 1 means mark as duplicate. The GBAM flags column has to be not compressed for patching to work as expected.
    /// For generating the markdup marks use: <time (../target/release/gbam_binary -v little.gbam > /tmp/testpipe.bam & samtools markdup -@5 -u /tmp/testpipe.bam /tmp/testoutpipe.bam & samtools view /tmp/testoutpipe.bam | awk '{print and($2, 0x400)!=0}' | ../target/release/gbam_binary --patch-gbam-with-dups little.gbam)>
    /// For pipes use <mkfifo> command.
    #[structopt(long)]
    patch_gbam_with_dups: bool,
    /// When sorting and converting file, only sort the indices of records but not the data itself.
    #[structopt(long)]
    index_sort: bool,
    /// Index file for use in Depth.
    #[structopt(long, parse(from_os_str))]
    index_file: Option<PathBuf>,
    /// Calculate uncompressed size of BAM file.
    #[structopt(long)]
    calc_uncompressed_size: bool,
}

/// Limited wrapper of `gbam_tools` converts BAM file to GBAM
/// file. Also limited tests may be run.
fn main() {
    let args = Cli::from_args();
    let arguments_strings: Vec<String> = env::args().collect();
    let full_command = arguments_strings.join(" ");
    if args.convert_to_gbam {
        convert(args, full_command);
    } else if args.test {
        test(args);
    } else if args.parallel_cigar_fetch {
        test_parallel_cigar_fetch(args);
    } else if args.depth {
        depth(args);
    } else if args.convert_to_bam {
        convert_to_bam(args);
    } else if args.flagstat {
        flagstat(args);
    } else if args.header {
        view_header(args);
    } else if args.view {
        let mut template = ParsingTemplate::new();
        template.set_all();
        view_file(args, template);
    } else if args.markdup_view {
        let mut template = ParsingTemplate::new();
        template.set_all_except(&[Fields::RawQual,Fields::RawSequence]);
        view_file(args, template);
    } else if args.patch_gbam_with_dups {
        patch_dups(args);
    }else if args.calc_uncompressed_size {
        test_file_uncompressed_size_fetch(args);
    }
}

fn convert(args: Cli, full_command: String) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path");
    let out_path = args
        .out_path
        .as_ref()
        .expect("Output path is mandatory for this operation.")
        .as_path()
        .to_str()
        .unwrap();
    if args.sort {
        bam_sort_to_gbam(in_path, out_path, Codecs::Lz4, args.sort_temp_mode, args.temp_dir, full_command, args.index_sort);
    } else {
        bam_to_gbam(in_path, out_path, Codecs::Lz4, full_command);
    }
}

fn convert_to_bam(args: Cli) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path.");
    let out_path = args
        .out_path
        .as_ref()
        .expect("Output path is mandatory for this operation.")
        .as_path()
        .to_str()
        .unwrap();
    gbam_to_bam(in_path, out_path);
}

fn flagstat(args: Cli) {
    let in_path = args
        .in_path
        .as_path()
        .to_str()
        .expect("Couldn't parse input path.");

    let file = File::open(in_path).unwrap();
    collect_stats(file);
}

fn test(args: Cli) {
    let mut tmplt = ParsingTemplate::new();
    tmplt.set(&Fields::RawCigar, true);

    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();

    let mut reader = Reader::new(file, tmplt).unwrap();
    let mut records = reader.records();
    let now = Instant::now();

    let mut u = 0;
    #[allow(unused_variables)]
    while let Some(rec) = records.next_rec() {
        u += base_coverage(&rec.cigar.as_ref().unwrap().0[..]);
    }
    println!("Record count {}", u);
    println!(
        "GBAM. Time elapsed querying POS and RAWCIGAR field throughout whole file: {}ms",
        now.elapsed().as_millis()
    );
    drop(records);
}

fn test_parallel_cigar_fetch(args: Cli) {
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();
    let temp_reader = Reader::new(file.try_clone().unwrap(), ParsingTemplate::new()).unwrap();
    let file_meta = temp_reader.file_meta;
    let total_records = temp_reader.amount;
    let now = Instant::now();
    
    (0..total_records).into_par_iter().chunks(500_000).for_each(|records_range| {
        let mut rec =  GbamRecord::default();
        let mut tmplt = ParsingTemplate::new();
        tmplt.set(&Fields::RawCigar, true);
    
        let mut reader = Reader::new_with_meta(file.try_clone().unwrap(), tmplt, &file_meta, None).unwrap();

        let mut collector = Vec::with_capacity(records_range.len());

        for rec_num in records_range {
            reader.fill_record(rec_num, &mut rec);
            collector.push(base_coverage(&rec.cigar.as_ref().unwrap().0[..]));
        }
    });

    println!(
        "Fetching CIGAR in parallel took: {}",
        now.elapsed().as_millis()
    );
}

fn test_file_uncompressed_size_fetch(args: Cli) {
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();

    let file_sz = file.metadata().unwrap().len();
    if file_sz == 0 {
        println!("File is empty.");
        return;
    }

    let mut reader = file;
    

    let mut buf: [u8; 1000] = [0; 1000];
    const OFFEST_IN_BGZF_FILE_TILL_BLOCK_SIZE_VALUE : usize = 128/8;
    let mut total_uncrompressed_size_of_file : usize = 0;
    const ERR : &str = "Couldn't parse the bgzf block.";
    loop {
        let cur_reader_pos = reader.seek(std::io::SeekFrom::Current(0)).unwrap();
        if file_sz == cur_reader_pos {
            break;
        }
        if file_sz-cur_reader_pos == 28 {
            break;
        }
        reader.read_exact(&mut buf[..OFFEST_IN_BGZF_FILE_TILL_BLOCK_SIZE_VALUE]).expect(ERR); 
        let block_size = reader.read_u16::<LittleEndian>().expect(ERR)+1;
        let uncompressed_info_start = cur_reader_pos+block_size as u64 - std::mem::size_of::<u32>() as u64;
        assert!(uncompressed_info_start < file_sz);
        reader.seek(std::io::SeekFrom::Start(uncompressed_info_start)).unwrap();
        let uncompressed_block_size = reader.read_u32::<LittleEndian>().expect(ERR);
        total_uncrompressed_size_of_file += uncompressed_block_size as usize;
        
    }

    println!("Total uncompressed size of file is: {}", total_uncrompressed_size_of_file);
}

fn read_index(index: PathBuf) -> Option<std::sync::Arc<Vec<u32>>> {
    let file = File::open(index).unwrap();
    let size = file.metadata().unwrap().len();
    let mut f = std::io::BufReader::new(file);
    let mut res = vec![0 as u32; (size/(std::mem::size_of::<u32>() as u64)).try_into().unwrap()];

    for slot in &mut res {
        *slot = f.read_u32::<LittleEndian>().unwrap();
    }

    Some(std::sync::Arc::new(res))
}

fn depth(args: Cli) {
    let in_path = args.in_path.as_path().to_str().unwrap();
    let gbam_file = File::open(in_path).unwrap();
    main_depth(gbam_file, args.bed_file.as_ref(), args.index_file.and_then(read_index), args.query, args.mapq, args.out_path, args.thread_num);
}

fn view_header(args: Cli){
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();
    let reader = Reader::new(file, ParsingTemplate::new()).unwrap();
    
    let header_len = (&reader.file_meta.get_sam_header()[..std::mem::size_of::<u32>()]).read_u32::<LittleEndian>().unwrap() as usize;
    let header_bytes = reader.file_meta.get_sam_header()[std::mem::size_of::<u32>()..std::mem::size_of::<u32>()+header_len].to_owned();
    let header = 
        String::from_utf8(header_bytes).unwrap();
   
    println!("{}", header);
}

fn parse_tag(tags: &[u8], target_tag: &str) -> Option<String> {
    // Convert the target tag to bytes
    let target_bytes = target_tag.as_bytes();
    let tag_length = target_bytes.len();

    if tag_length == 0 || tag_length + 1 >= tags.len() {
        // Ensure the target tag is not empty and fits within the input slice
        return None;
    }

    let mut i = 0;
    while i <= tags.len() - tag_length - 1 {
        // Look for the target tag followed by its type identifier
        if tags[i..i + tag_length] == target_bytes[..] {
            // Skip the tag name and type
            i += tag_length + 1;
            let mut result = Vec::new();

            // Read until null terminator or end of tags
            while i < tags.len() && tags[i] != 0 {
                result.push(tags[i]);
                i += 1;
            }

            if !result.is_empty() {
                return String::from_utf8(result).ok();
            }
        }
        i += 1;
    }
    None
}


fn view_file(args: Cli, template: ParsingTemplate){
    let file = File::open(args.in_path.as_path().to_str().unwrap()).unwrap();

    let mut reader = Reader::new_with_index(file, template, args.index_file.and_then(read_index)).unwrap();

    let st = std::io::stdout();
    let lock = st.lock();
    let mut stdout = BufWriter::with_capacity(64 * 1024, lock);

    const BAM_MAGIC: &[u8; 4] = b"BAM\x01";
    stdout.write_all(BAM_MAGIC).unwrap();
    stdout.write_all(reader.file_meta.get_sam_header()).unwrap();
    
    let mut records = reader.records();
    let mut buf = Vec::new();
    while let Some(rec) = records.next_rec() {
        rec.convert_to_bytes(&mut buf);
        if stdout.write_all(&buf).is_err() {
            break;
        }
    }
}



fn patch_dups(args: Cli){

    let file = OpenOptions::new()
        .write(true)
        .read(true)
        .open(args.in_path.as_path().to_str().unwrap())
        .unwrap();

    let reader = Reader::new_with_index(file.try_clone().unwrap(), ParsingTemplate::new(), args.index_file.and_then(read_index)).unwrap();
    let file_meta = reader.file_meta.clone();

    let mut buf = Vec::new();

    let codec = file_meta.get_field_codec(&Fields::Flags);
    assert!(codec == &Codecs::NoCompression);

    let mut read_manual = BufReader::with_capacity(MEGA_BYTE_SIZE, file.try_clone().unwrap());
    let mut write_manual = BufWriter::with_capacity(MEGA_BYTE_SIZE, file.try_clone().unwrap());
    for block in file_meta.view_blocks(&Fields::Flags){
        let available_in_block = block.numitems;
        buf.resize(block.block_size as usize, 0);
        read_manual.seek(SeekFrom::Start(block.seekpos)).unwrap();
        read_manual.read_exact(&mut buf).unwrap();
        let slice = &mut buf[..];
        for (chunk, is_dup) in zip_eq(slice.chunks_mut(2), std::io::stdin().lock().lines().take(available_in_block as usize)){
            let mut val = (&chunk[..]).read_u16::<byteorder::LittleEndian>().unwrap();
            if is_dup.unwrap() == "1" {
                val = val | 0x400;
            }
            (&mut chunk[..]).write_u16::<byteorder::LittleEndian>(val).unwrap();
        }
        write_manual.seek(SeekFrom::Start(block.seekpos)).unwrap();
        write_manual.write_all(&buf).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag_with_exact_match() {
        let tags = b"TAGZvalue1\0LONGERTAGZvalue2\0";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, Some("value1".to_string()));
    }

    #[test]
    fn test_parse_tag_with_longer_tag() {
        let tags = b"TAGZvalue1\0LONGERTAGZvalue2\0";
        let result = parse_tag(tags, "LONGERTAG");
        assert_eq!(result, Some("value2".to_string()));
    }

    #[test]
    fn test_parse_tag_with_nonexistent_tag() {
        let tags = b"TAGZvalue1\0LONGERTAGZvalue2\0";
        let result = parse_tag(tags, "NONEXISTENT");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_tag_with_empty_tags() {
        let tags = b"";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_tag_with_empty_target_tag() {
        let tags = b"TAGZvalue1\0LONGERTAGZvalue2\0";
        let result = parse_tag(tags, "");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_tag_with_multiple_tags() {
        let tags = b"TAGZvalue1\0TAGZvalue2\0";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, Some("value1".to_string()));
    }

    #[test]
    fn test_parse_tag_with_null_terminator_at_end() {
        let tags = b"TAGZvalue1\0";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, Some("value1".to_string()));
    }

    #[test]
    fn test_parse_tag_with_partial_tag_match() {
        let tags = b"TAZvalue1\0";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_tag_with_no_null_terminator() {
        let tags = b"TAGZvalue1";
        let result = parse_tag(tags, "TAG");
        assert_eq!(result, Some("value1".to_string()));
    }
}
