from enum import Enum, auto
from typing import List, Tuple, Optional, BinaryIO
import pysam


class FieldType(Enum):
    FIXED_SIZED = auto()
    VARIABLE_SIZED = auto()


class Fields(Enum):
    REF_ID = auto()
    POS = auto()
    MAPQ = auto()
    BIN = auto()
    FLAGS = auto()
    NEXT_REF_ID = auto()
    NEXT_POS = auto()
    TEMPLATE_LENGTH = auto()
    READ_NAME = auto()
    RAW_CIGAR = auto()
    RAW_SEQUENCE = auto()
    RAW_QUAL = auto()
    RAW_TAGS = auto()
    LNAME = auto()
    NCIGAR = auto()
    SEQUENCE_LENGTH = auto()
    RAW_TAGS_LEN = auto()
    RAW_SEQ_LEN = auto()

    @staticmethod
    def iterator():
        return iter(Fields)

    @staticmethod
    def is_data_field(field):
        return field in {
            Fields.REF_ID, Fields.POS, Fields.MAPQ, Fields.BIN, Fields.FLAGS, Fields.NEXT_REF_ID,
            Fields.NEXT_POS, Fields.TEMPLATE_LENGTH, Fields.READ_NAME, Fields.RAW_CIGAR,
            Fields.RAW_SEQUENCE, Fields.RAW_QUAL, Fields.RAW_TAGS
        }

    @staticmethod
    def field_type(field):
        fixed_sized_fields = {
            Fields.REF_ID, Fields.POS, Fields.LNAME, Fields.MAPQ, Fields.BIN, Fields.NCIGAR,
            Fields.FLAGS, Fields.SEQUENCE_LENGTH, Fields.NEXT_REF_ID, Fields.NEXT_POS,
            Fields.TEMPLATE_LENGTH, Fields.RAW_SEQ_LEN, Fields.RAW_TAGS_LEN
        }
        return FieldType.FIXED_SIZED if field in fixed_sized_fields else FieldType.VARIABLE_SIZED


class Column:
    """Base class for all columns."""
    
    def __init__(self, field: Fields, stat_collector: Optional[object] = None):
        self.field = field
        self.stat_collector = stat_collector
        self.data = []  # Placeholder for storing column data

    def add_data(self, value):
        """Adds data to the column."""
        self.data.append(value)

    def compress(self):
        """Placeholder for compression logic."""
        raise NotImplementedError("Compression should be implemented in subclasses")


class FixedColumn(Column):
    """Column for fixed-size fields."""
    
    def __init__(self, field: Fields, stat_collector: Optional[object] = None):
        super().__init__(field, stat_collector)
    
    def compress(self):
        """Compression logic for fixed-size data."""
        return b"".join(value.to_bytes(4, 'little') if isinstance(value, int) else value for value in self.data)


class VariableColumn(Column):
    """Column for variable-sized fields."""
    
    def __init__(self, field: Fields, stat_collector: Optional[object] = None):
        super().__init__(field, stat_collector)
    
    def compress(self):
        """Compression logic for variable-sized data."""
        return b"".join(len(value).to_bytes(4, 'little') + value.encode() if isinstance(value, str) else value for value in self.data)


class Compressor:
    """Handles multi-threaded compression."""
    
    def __init__(self, thread_num: int):
        self.thread_num = thread_num

    def compress_column(self, column: Column):
        """Compresses a single column (placeholder logic)."""
        return column.compress()


class FileMeta:
    """Metadata for the file, including codecs and references."""
    
    def __init__(self, codec, ref_seqs, sam_header):
        self.codec = codec
        self.ref_seqs = ref_seqs
        self.sam_header = sam_header


class FileInfo:
    """Holds file information including version, column counts, and sorting status."""
    
    def __init__(self, version, field_count, column_count, full_command, is_sorted):
        self.version = version
        self.field_count = field_count
        self.column_count = column_count
        self.full_command = full_command
        self.is_sorted = is_sorted


class Writer:
    """Writer class for handling BAM-like structured files with compression."""
    
    def __init__(
        self,
        inner: BinaryIO,
        codecs: List[str],  
        thread_num: int,
        collect_stats_for: List[Fields],
        ref_seqs: List[Tuple[str, int]],
        sam_header: bytes,
        full_command: str,
        is_sorted: bool,
    ):
        inner.seek(64)  # FILE_INFO_SIZE equivalent

        self.columns = []
        count = 0

        for field in Fields.iterator():
            if Fields.is_data_field(field):
                stat_collector = {} if field in collect_stats_for else None
                
                if Fields.field_type(field) == FieldType.FIXED_SIZED:
                    col = FixedColumn(field, stat_collector)
                else:
                    count += 1  # Index column adjustment
                    col = VariableColumn(field, stat_collector)

                self.columns.append(col)
                count += 1

        assert count == len(list(Fields.iterator()))  # Equivalent to `debug_assert!(count == FIELDS_NUM);`

        self.file_meta = FileMeta(codecs[0], ref_seqs, sam_header)
        self.inner = inner
        self.compressor = Compressor(thread_num)
        self.file_info = FileInfo([1, 0], 0, 0, full_command, is_sorted)

def extract_fields(record) -> dict:
    """Extracts relevant fields from a BAM record."""
    return {
        Fields.REF_ID: record.reference_id,
        Fields.POS: record.reference_start,
        Fields.MAPQ: record.mapping_quality,
        Fields.BIN: record.bin,
        Fields.FLAGS: record.flag,
        Fields.NEXT_REF_ID: record.next_reference_id,
        Fields.NEXT_POS: record.next_reference_start,
        Fields.TEMPLATE_LENGTH: record.template_length,
        Fields.READ_NAME: record.query_name,
        Fields.RAW_CIGAR: record.cigarstring if record.cigarstring else "",
        Fields.RAW_SEQUENCE: record.query_sequence if record.query_sequence else "",
        Fields.RAW_QUAL: record.query_qualities.tobytes() if record.query_qualities else b"",
        Fields.RAW_TAGS: record.get_tags(),
    }


def convert_bam_to_gbam(bam_file: str, gbam_file: str):
    """Reads a BAM file and writes it to GBAM format using the Writer class."""
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(gbam_file, "wb") as gbam:
        
        ref_seqs = [(ref, length) for ref, length in zip(bam.references, bam.lengths)]
        sam_header = bam.header.to_dict()
        codecs = ["zlib"]  # Placeholder codec, can be changed
        thread_num = 4  # Example thread count
        collect_stats_for = [Fields.REF_ID, Fields.POS]
        full_command = f"convert_bam_to_gbam {bam_file} {gbam_file}"
        is_sorted = bam.header.get("HD", {}).get("SO", "") == "coordinate"
        
        writer = Writer(
            gbam,
            codecs,
            thread_num,
            collect_stats_for,
            ref_seqs,
            str(sam_header).encode(),
            full_command,
            is_sorted,
        )
        
        for record in bam.fetch():
            extracted_data = extract_fields(record)
            for column in writer.columns:
                value = extracted_data.get(column.field)
                if value is not None:
                    column.add_data(value)
        
        for column in writer.columns:
            compressed_data = writer.compressor.compress_column(column)
            gbam.write(compressed_data)


if __name__ == "__main__":
    input_bam = "input.bam"
    output_gbam = "output.gbam"
    convert_bam_to_gbam(input_bam, output_gbam)