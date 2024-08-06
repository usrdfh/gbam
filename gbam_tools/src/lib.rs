//! This crate contains tools for interacting with GBAM file format.
#![allow(missing_docs)]

use std::mem;

pub mod bam {
    /// BAM to GBAM converter
    pub mod bam_to_gbam;
    /// GBAM to BAM converter
    pub mod gbam_to_bam;
}
///
pub mod utils {
    /// BED reader
    pub mod bed;
}

pub mod reader {
    pub mod column;
    pub mod parse_tmplt;
    /// GBAM reader
    #[allow(clippy::module_inception)]
    pub mod reader;
    pub mod record;
    pub mod records;

}

#[cfg(not(feature = "python-ffi"))]
pub mod query {
    pub mod cigar;
    pub mod depth;
    pub mod flagstat;
    pub mod int2str;
    pub mod markdup {
        pub mod markdup;
        mod sorted_storage;
    }
}


/// Manages parallel compression
mod compressor;
/// Meta information for GBAM file
pub mod meta;
/// Manages stats collection
mod stats;
/// GBAM writer
pub mod writer;

// use self::writer::Writer;
// pub use {ParsingTemplate, Reader};
use self::writer::Writer;
pub use bam::bam_to_gbam::{bam_sort_to_gbam, bam_to_gbam};
pub use meta::Codecs;
pub use bam_tools::record::fields::Fields;


const U32_SIZE: usize = mem::size_of::<u32>();
const MEGA_BYTE_SIZE: usize = 1_048_576;

/// 16777216 bytes
const SIZE_LIMIT: usize = 8 * MEGA_BYTE_SIZE;
static GBAM_MAGIC: &[u8] = b"geeBAM10";
