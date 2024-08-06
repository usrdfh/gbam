use crate::{U16_SIZE, U32_SIZE, U8_SIZE};
/// This module contains declaration of Fields enum, which is used to ease
/// manipulating BAM record fields. Some of the fields are not in the BAM spec
/// and used exclusively in GBAM crate.
use serde::{Deserialize, Serialize};

use std::slice::Iter;

// To avoid visual clutter (no need to write Fields::* each time).
use self::Fields::*;

pub const FIELDS_NUM: usize = 18;
/// Fields which contain data (not index fields).
#[allow(dead_code)]
pub const DATA_FIELDS_NUM: usize = 13;
/// Types of fields contained in BAM file.
/// WARNING: ORDER IS IMPORTANT! 
/// gbam_tools/src/reader/reader.rs:85, if the order is incorrect, bam1_t structure data field
/// will be filled in a wrong order and it would break it.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[allow(missing_docs)]
pub enum Fields {
    /// Data fields
    #[allow(clippy::upper_case_acronyms)]
    RefID,
    Pos,
    Mapq,
    Bin,
    Flags,
    #[allow(clippy::upper_case_acronyms)]
    NextRefID,
    NextPos,
    TemplateLength,
    ReadName,
    RawCigar,
    RawSequence,
    RawQual,
    RawTags,
    /// Index fields
    LName,
    NCigar,
    SequenceLength,
    RawTagsLen, // Not in BAM spec, needed for index GBAM file
    RawSeqLen,  // Not in BAM spec, needed for index GBAM file
}

impl Fields {
    /// Returns iterator over enum fields
    pub fn iterator() -> Iter<'static, Fields> {
        static FIELDS: [Fields; FIELDS_NUM] = [
            // Data fields
            RefID,
            Pos,
            Mapq,
            Bin,
            Flags,
            NextRefID,
            NextPos,
            TemplateLength,
            ReadName,
            RawCigar,
            RawSequence,
            RawQual,
            RawTags,
            // Index fields
            LName,
            NCigar,
            SequenceLength,
            RawSeqLen,
            RawTagsLen,
        ];
        FIELDS.iter()
    }
}

/// Fields holding index are not data fields
pub fn is_data_field(field: &Fields) -> bool {
    matches!(
        field,
        RefID
            | Pos
            | Mapq
            | Bin
            | Flags
            | NextRefID
            | NextPos
            | TemplateLength
            | ReadName
            | RawCigar
            | RawSequence
            | RawQual
            | RawTags
    )
}

/// Variable sized fields have fixed size value NONE
/// This is **not** BAM sizes!
pub fn field_item_size(field: &Fields) -> Option<usize> {
    match field {
        RefID => Some(U32_SIZE),
        Pos => Some(U32_SIZE),
        LName => Some(U32_SIZE),
        Mapq => Some(U8_SIZE),
        Bin => Some(U16_SIZE),
        NCigar => Some(U32_SIZE),
        Flags => Some(U16_SIZE),
        SequenceLength => Some(U32_SIZE),
        NextRefID => Some(U32_SIZE),
        NextPos => Some(U32_SIZE),
        TemplateLength => Some(U32_SIZE),
        RawTagsLen => Some(U32_SIZE),
        RawSeqLen => Some(U32_SIZE),
        ReadName => None,
        RawCigar => None,
        RawSequence => None,
        RawQual => None,
        RawTags => None,
    }
}

impl std::fmt::Display for Fields {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}

/// Type of Field.
#[derive(Debug)]
#[allow(missing_docs)]
pub enum FieldType {
    VariableSized,
    FixedSized,
}

// May be useful:
pub fn field_type(field: &Fields) -> FieldType {
    match field {
        // Fixed size fields
        Fields::RefID
        | Fields::Pos
        | Fields::LName
        | Fields::Mapq
        | Fields::Bin
        | Fields::NCigar
        | Fields::Flags
        | Fields::SequenceLength
        | Fields::NextRefID
        | Fields::NextPos
        | Fields::TemplateLength
        | Fields::RawSeqLen
        | Fields::RawTagsLen => FieldType::FixedSized,
        // Variable size fields
        Fields::ReadName
        | Fields::RawCigar
        | Fields::RawSequence
        | Fields::RawQual
        | Fields::RawTags => FieldType::VariableSized,
    }
}

/// Returns enum name of index field for particular variable sized field
pub fn var_size_field_to_index(field: &Fields) -> Fields {
    match field {
        Fields::ReadName => Fields::LName,
        Fields::RawQual => Fields::SequenceLength,
        Fields::RawSequence => Fields::RawSeqLen,
        Fields::RawTags => Fields::RawTagsLen,
        Fields::RawCigar => Fields::NCigar,
        _ => panic!("Unreachable"),
    }
}
