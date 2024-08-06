use gbam_tools::reader::{parse_tmplt::ParsingTemplate, reader::Reader, record::GbamRecord};
use gbam_tools::Fields;
use libc::{uint64_t, uint8_t};
use rust_htslib::htslib::{bam1_t, kstring_t, sam_hdr_destroy, sam_hdr_parse, sam_hdr_t, uint64_u};
use byteorder::{LittleEndian, ReadBytesExt};
use std::ffi::CStr;
use std::fs::File;
use std::mem::MaybeUninit;
use std::path::Path;

const DISABLE_READNAME : u64 = 1<<0;
const DISABLE_CIGAR : u64 = 1<<1;
const DISABLE_SEQ : u64 = 1<<2;
const DISABLE_QUAL : u64 = 1<<3;
const DISABLE_TAGS : u64 = 1<<4;

#[no_mangle]
pub extern "C" fn get_gbam_reader(gbam_path: *const libc::c_char, flags: u64) -> *mut Reader {
    let gbam_path_str: String = unsafe {
        assert!(!gbam_path.is_null());
        CStr::from_ptr(gbam_path).to_string_lossy().into_owned()
    };



    let gbam_path = Path::new(&gbam_path_str);

    if !gbam_path.exists() {
        eprintln!(
            "The file doesnt exists or couldn't be found by Rust code: {}",
            gbam_path.to_str().unwrap()
        );
        return std::ptr::null_mut();
    }

    let file = File::open(gbam_path).unwrap();
    let mut tmplt = ParsingTemplate::new();
    let mapping = vec![
        (DISABLE_READNAME, Fields::ReadName),
        (DISABLE_CIGAR, Fields::RawCigar),
        (DISABLE_SEQ, Fields::RawSequence),
        (DISABLE_QUAL, Fields::RawQual),
        (DISABLE_TAGS, Fields::RawTags),
    ];

    let mut field_to_not_set = Vec::new();
    for (f, field) in mapping {
        if (f&flags) > 0 {
            field_to_not_set.push(field);
        }
    }

    tmplt.set_all_except(&field_to_not_set[..]);
    let reader = Reader::new(file, tmplt).unwrap();

    return Box::into_raw(Box::new(reader));
}

#[no_mangle]
pub extern "C" fn get_records_num(reader: *mut Reader) -> u64 {
    unsafe {
        (*reader).amount as u64
    }
}

#[no_mangle]
pub extern "C" fn get_bam_record(reader: *mut Reader, rec_num: u64) -> bam1_t {
    let mut rec = GbamRecord::default();
    unsafe {
        (*reader).fill_record(rec_num as usize, &mut rec);
    }

    return rec.get_hts_repr();
}



// Returns 0 on success,
//        -1 on EOF,
#[no_mangle]
pub extern "C" fn fill_next_record(reader: *mut Reader, rec: *mut bam1_t) -> i32 {
    unsafe {
        if (*reader).amount == (*reader).next_not_read {
            return -1;
        }
        let mut rec = Box::from_raw(rec);
        
        (*reader).fill_bam1_t_record((*reader).next_not_read, &mut rec);
        (*reader).next_not_read += 1;

        Box::into_raw(rec);
    }

    return 0;
}

#[no_mangle]
pub extern "C" fn get_empty_bam1_t() -> *mut bam1_t {
    let mut v : bam1_t = unsafe { MaybeUninit::uninit().assume_init() };
    
    v.core.l_qname = 0;
    v.core.l_qseq = 0;
    v.core.n_cigar = 0;

    let mut inner_vec = Vec::<u8>::new();
    inner_vec.resize(1, b'\0');


    let capacity = inner_vec.capacity();
    let len = inner_vec.len();
    let data = inner_vec.into_boxed_slice();
    unsafe {
        let ptr = Box::into_raw(data);
        v.data = ptr.as_mut().unwrap().as_mut_ptr();
        v.l_data = len as i32;
        v.m_data = capacity as u32;
    }

    Box::into_raw(Box::new(v))
}


#[no_mangle]
pub extern "C" fn get_bam_record_fast(reader: *mut Reader, rec_num: u64) -> *mut bam1_t {
    let ptr = get_empty_bam1_t();
    unsafe {

        let mut rec = Box::from_raw(ptr);
       

        (*reader).fill_bam1_t_record(rec_num as usize, &mut rec);
        
        Box::into_raw(rec);
    }
    ptr
}

#[no_mangle]
pub extern "C" fn free_bam_record(rec: *mut bam1_t) {
    unsafe {
        let rec = Box::from_raw(rec);
        drop(Box::from_raw(rec.data));
        drop(rec);
    }
}

#[no_mangle]
pub extern "C" fn free_gbam_reader(reader: *mut Reader) {
    unsafe {
        drop(Box::<Reader>::from_raw(reader));
    }
}

#[no_mangle]
pub extern "C" fn get_header(reader: *mut Reader) -> *mut sam_hdr_t {
    unsafe {
        let header_len = (&(*reader).file_meta.get_sam_header()[..std::mem::size_of::<u32>()]).read_u32::<LittleEndian>().unwrap() as usize;
        let header_bytes = (*reader).file_meta.get_sam_header()[std::mem::size_of::<u32>()..std::mem::size_of::<u32>()+header_len].to_owned();

        sam_hdr_parse(header_len as u64, header_bytes.as_ptr() as *const i8)
    }
}


#[no_mangle]
pub extern "C" fn free_header(hdr: *mut sam_hdr_t) {
    unsafe{
        sam_hdr_destroy(hdr);
    }   
}

// int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b) HTS_RESULT_USED;
// Make a python test
// Need to handle the header
// Verify in Python that the C API works fine and there are no leaks
// Integrate this into SAMTOOLS view code and verify it works
// Integrate this into SAMTOOLS depth and SAMTOOLS markdup


// The serializer from bincode does bullshit, it's just for saving not for exchangin data with samtools
// Pass the data to samtools directly
// The problem to figure out is the header, it has to be encoded and passed to samtools somehow too


// Create a header from existing text.
//
// @param l_text    Length of text
// @param text      Header text
// @return A populated sam_hdr_t structure on success; NULL on failure.
// @note The text field of the returned header will be NULL, and the l_text
// field will be zero.
//
// The sam_hdr_t struct returned by a successful call should be freed
// via sam_hdr_destroy() when it is no longer needed.
//
// HTSLIB_EXPORT
// sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text);
// use rust_htslib::htslib::{bam_write1, sam_hdr_t, sam_hdr_parse};
