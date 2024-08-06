use super::{fields::Fields, tags::get_tag};
use crate::{U16_SIZE, U32_SIZE, U8_SIZE};

use byteorder::{ByteOrder, LittleEndian, ReadBytesExt, WriteBytesExt};
use std::borrow::Cow;
use std::ops::{Deref, DerefMut};
use std::io::Write;

use super::tags::get_hit_count;

/// Provides convenient access to BAM-style raw read (record bytes)
/// Cow is used so BAMRawRecord can either own or borrow underlying data (if it won't be mutated).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BAMRawRecord<'a>(pub Cow<'a, [u8]>);

impl<'a> BAMRawRecord<'a> {
    /// Change size of record
    pub fn resize(&mut self, new_len: usize) {
        self.0.to_mut().resize(new_len, Default::default());
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
            Fields::RawCigar => self.get_cigar(get_cigar_offset()),
            Fields::RawSequence => self.get_slice(get_seq_offset(), self.get_var_field_len(field)),
            Fields::RawQual => self.get_slice(get_qual_offset(), self.l_seq() as usize),
            Fields::RawTags => self.get_slice(get_tags_offset(), self.0.len() - get_tags_offset()),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    pub(crate) fn get_name(&self) -> &str {
        unsafe { std::str::from_utf8_unchecked(self.get_bytes(&Fields::ReadName)) }
    }

    pub(crate) fn get_refid(&self) -> i32 {
        self.get_bytes(&Fields::RefID)
            .read_i32::<LittleEndian>()
            .unwrap()
    }

    // Calculates range of bytes containing specified field.
    pub fn get_range(&self, field: &Fields) -> std::ops::Range<usize> {
        match field {
            Fields::ReadName => 32..(32 + self.get_var_field_len(field)),
            _ => panic!("This field is not supported: {} \n", *field as usize),
        }
    }

    /// Extracts CIGAR from tags if it didn't fit into CIGAR field
    fn get_cigar(&self, cigar_offset: usize) -> &[u8] {
        let ref_id = self
            .get_bytes(&Fields::RefID)
            .read_i32::<LittleEndian>()
            .unwrap();
        let pos = self
            .get_bytes(&Fields::Pos)
            .read_i32::<LittleEndian>()
            .unwrap();
        if ref_id < 0 || pos < 0 || self.get_var_field_len(&Fields::RawCigar) == 0 {
            return &[];
        }
        let cigar_field_data =
            self.get_slice(cigar_offset, self.get_var_field_len(&Fields::RawCigar));

        let mut first_op_bytes = &cigar_field_data[..U32_SIZE];
        let first_op = first_op_bytes.read_u32::<LittleEndian>().unwrap() as usize;
        let mut n_cigar_bytes = self.get_bytes(&Fields::NCigar);
        let n_cigar = n_cigar_bytes.read_u16::<LittleEndian>().unwrap() as usize;

        if (first_op & 0xf) != 4
            || (first_op >> 4) != self.get_var_field_len(&Fields::RawSequence)
            || n_cigar != 2
        {
            return cigar_field_data;
        }
        let cigar_tag = &[b'C', b'G'];
        match self.get_tag(cigar_tag) {
            Some(cigar) => cigar,
            None => cigar_field_data, //panic!("CIGAR in tags not found!"),
        }
    }

    pub(crate) fn get_tag(&self, tag: &[u8; 2]) -> Option<&[u8]> {
        get_tag(self.get_bytes(&Fields::RawTags), tag).map(|tag_val| tag_val.0)
    }

    pub fn get_hit_count(&self) -> Option<i32> {
        get_hit_count(self.get_bytes(&Fields::RawTags))
    }
}

impl<'a> From<Vec<u8>> for BAMRawRecord<'a> {
    fn from(bytes: Vec<u8>) -> Self {
        Self(Cow::Owned(bytes))
    }
}

// Source: https://github.com/zaeleus/noodles/blob/316ec6f42960e4540bb2acc45b5653fb00b9970c/noodles-bam/src/record.rs#L324
impl<'a> Default for BAMRawRecord<'a> {
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

impl<'a> Deref for BAMRawRecord<'a> {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl<'a> DerefMut for BAMRawRecord<'a> {
    fn deref_mut(&mut self) -> &mut [u8] {
        self.0.to_mut()
    }
}

fn get_cig_op(val: u32) -> char {
    match val {
        0 => 'M',
        1 => 'I',
        2 => 'D',
        3 => 'N',
        4 => 'S',
        5 => 'H',
        6 => 'P',
        7 => '=',
        8 => 'X',
        _ => panic!("Value <{}> is not supported.", val),
    }
}

/// Decodes CIGAR bytes into a string
pub fn decode_cigar(bytes: &[u8]) -> String {
    let mut res = String::new();
    res.reserve(3 * bytes.len());
    for cig in bytes
        .chunks(U32_SIZE)
        .map(|mut slice| slice.read_u32::<LittleEndian>().unwrap())
    {
        let op = cig & 0xf;
        let op_len = cig >> 4;
        res.push_str(&op_len.to_string());
        res.push(get_cig_op(op));
    }
    res
}

fn get_seq_base(val: u32) -> char {
    match val {
        0 => '=',
        1 => 'A',
        2 => 'C',
        3 => 'M',
        4 => 'G',
        5 => 'R',
        6 => 'S',
        7 => 'V',
        8 => 'T',
        9 => 'W',
        10 => 'Y',
        11 => 'H',
        12 => 'K',
        13 => 'D',
        14 => 'B',
        15 => 'N',
        _ => panic!("Value <{}> is not supported.", val),
    }
}

/// Decodes sequence bytes into string
pub fn decode_seq(bytes: &[u8], res: &mut String) {
    res.clear();
    for byte in bytes {
        let first = (byte >> 4) as u32;
        let second = (byte & 0xf) as u32;
        res.push(get_seq_base(first));
        // WARNING TODO: Assumption that the last half will contain value or zero (it is not mandated.)
        if second != 0 {
            res.push(get_seq_base(second));
        }
    }
}


/// From NOODLES crate.
pub fn put_sequence<W: Write>(mut dst: &mut W, read_length: usize, sequence: &String) -> std::io::Result<()>

{
    if sequence.is_empty() {
        return Ok(());
    }

    // § 1.4.10 "`SEQ`" (2022-08-22): "If not a '*', the length of the sequence must equal the sum
    // of lengths of `M`/`I`/`S`/`=`/`X` operations in `CIGAR`."
    if read_length > 0 && sequence.len() != read_length {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "read length-sequence length mismatch",
        ));
    }

    let mut bases = sequence.chars();

    while let Some(l) = bases.next() {
        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these as zero."
        let r = bases.next().unwrap_or('\0');
        let b = encode_base(l) << 4 | encode_base(r);
        dst.write_u8(b)?;
    }

    Ok(())
}

fn encode_base(base: char) -> u8 {
    match base {
        '=' => 0,
        'A' => 1,
        'C' => 2,
        'M' => 3,
        'G' => 4,
        'R' => 5,
        'S' => 6,
        'V' => 7,
        'T' => 8,
        'W' => 9,
        'Y' => 10,
        'H' => 11,
        'K' => 12,
        'D' => 13,
        'B' => 14,
        // § 4.2.3 SEQ and QUAL encoding (2021-06-03): "The case-insensitive base codes ... are
        // mapped to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
        _ => 15,
    }
}
