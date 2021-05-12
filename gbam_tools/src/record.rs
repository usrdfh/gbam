use super::{Fields, U16_SIZE, U32_SIZE, U8_SIZE};
use byteorder::{ByteOrder, LittleEndian};
use std::ops::{Deref, DerefMut};

/// Provides convenient access to record bytes
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct RawRecord(pub Vec<u8>);

impl RawRecord {
    /// Change size of record
    pub fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    fn get_slice(&self, offset: usize, len: usize) -> &[u8] {
        &self.0[offset..offset + len]
    }

    fn l_read_name(&self) -> u8 {
        self.get_bytes(&Fields::LName)[0]
    }

    fn n_cigar_op(&self) -> u16 {
        LittleEndian::read_u16(self.get_bytes(&Fields::NCigar))
    }

    fn l_seq(&self) -> u32 {
        LittleEndian::read_u32(self.get_bytes(&Fields::SequenceLength))
    }

    /// Values of fields containg length of other fields
    pub fn get_len_val(&self, field: &Fields) -> usize {
        match field {
            Fields::LName => self.l_read_name() as usize,
            Fields::SequenceLength => self.l_seq() as usize,
            Fields::NCigar => self.n_cigar_op() as usize,
            Fields::RawTagsLen => self.0.len() - self.get_offset(&Fields::RawTags),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    /// Calculates actual size of variable length field in bytes.
    pub fn get_var_field_len(&self, field: &Fields) -> usize {
        match field {
            Fields::ReadName => self.l_read_name() as usize,
            Fields::RawCigar => U32_SIZE * self.n_cigar_op() as usize,
            Fields::RawSequence => ((self.l_seq() + 1) / 2) as usize,
            Fields::RawQual => self.l_seq() as usize,
            Fields::RawTags => self.0.len() - self.get_offset(&Fields::RawTags),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    fn get_offset(&self, field: &Fields) -> usize {
        match field {
            Fields::ReadName => 32,
            Fields::RawCigar => {
                self.get_offset(&Fields::ReadName) + self.get_var_field_len(&Fields::ReadName)
            }
            Fields::RawSequence => {
                self.get_offset(&Fields::RawCigar) + self.get_var_field_len(&Fields::RawCigar)
            }
            Fields::RawQual => {
                self.get_offset(&Fields::RawSequence) + self.get_var_field_len(&Fields::RawSequence)
            }
            Fields::RawTags => {
                self.get_offset(&Fields::RawQual) + self.get_var_field_len(&Fields::RawQual)
            }
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    /// Returns bytes of specified field
    pub fn get_bytes(&self, field: &Fields) -> &[u8] {
        let get_cigar_offset = || -> usize { (32 + self.l_read_name()) as usize };
        let get_seq_offset =
            || -> usize { get_cigar_offset() + U32_SIZE * self.n_cigar_op() as usize };
        let get_qual_offset = || -> usize { get_seq_offset() + ((self.l_seq() + 1) / 2) as usize };
        let get_tags_offset = || -> usize { get_qual_offset() + self.l_seq() as usize };
        match field {
            Fields::RefID => self.get_slice(0, U32_SIZE),
            Fields::Pos => self.get_slice(4, U32_SIZE),
            Fields::LName => self.get_slice(8, U8_SIZE),
            Fields::Mapq => self.get_slice(9, U8_SIZE),
            Fields::Bin => self.get_slice(10, U16_SIZE),
            Fields::NCigar => self.get_slice(12, U16_SIZE),
            Fields::Flags => self.get_slice(14, U16_SIZE),
            Fields::SequenceLength => self.get_slice(16, U32_SIZE),
            Fields::NextRefID => self.get_slice(20, U32_SIZE),
            Fields::NextPos => self.get_slice(24, U32_SIZE),
            Fields::TemplateLength => self.get_slice(28, U32_SIZE),
            Fields::ReadName => self.get_slice(32, self.get_var_field_len(field)),
            Fields::RawCigar => self.get_slice(get_cigar_offset(), self.get_var_field_len(field)),
            Fields::RawSequence => self.get_slice(get_seq_offset(), self.get_var_field_len(field)),
            Fields::RawQual => self.get_slice(get_qual_offset(), self.l_seq() as usize),
            Fields::RawTags => self.get_slice(get_tags_offset(), self.0.len() - get_tags_offset()),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }
}

impl From<Vec<u8>> for RawRecord {
    fn from(bytes: Vec<u8>) -> Self {
        Self(bytes)
    }
}

// Source: https://github.com/zaeleus/noodles/blob/316ec6f42960e4540bb2acc45b5653fb00b9970c/noodles-bam/src/record.rs#L324
impl Default for RawRecord {
    fn default() -> Self {
        Self::from(vec![
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ])
    }
}

impl Deref for RawRecord {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl DerefMut for RawRecord {
    fn deref_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}
