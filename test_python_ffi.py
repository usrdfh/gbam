import gbam_tools
from tempfile import NamedTemporaryFile
import pysam
from enum import IntEnum


# These values should have same values as in Rust version
class Fields(IntEnum):
    REFID = 0
    POS = 1
    MAPQ = 2 
    BIN = 3
    FLAGS = 4
    NEXTREFID = 5
    NEXTPOS = 6
    READNAME = 7
    RAWCIGAR = 8
    RAWSEQUENCE = 9
    RAWQUAL = 10
    RAWTAGS = 11



FIELDS_NUM = 12

def convert(bam_path, gbam_path):
    from gbam_tools import bam_to_gbam_python
    bam_to_gbam_python(bam_path, gbam_path)


def get_reader(path, parsing_tmplt):
    from gbam_tools import Reader, ParsingTemplate
    return Reader(path, parsing_tmplt)


# Converts BAM file to GBAM file, performs tests, deletes GBAM file.
def test(args):
    # input_path = args.input_path

    # output_file = NamedTemporaryFile()
    # convert(input_path, output_file.name)

    # test_combinations = 0

    # for field_to_omit in list(map(int, Fields)):
    print("Here")
    bam_path = "../test_data/1.bam"
    gbam_path = "../test_data/res.gbam"
    check_if_equal(bam_path, gbam_path)

    print("TESTS PASSED.")

def get_parsing_tmpl(fields_to_parse):
    from gbam_tools import ParsingTemplate
    fields = [False] * FIELDS_NUM
    for field in fields_to_parse:
        fields[field] = True
    return ParsingTemplate(fields)

# Checks if the data in both BAM and GBAM files is equal.


def check_if_equal(bam_path, gbam_path, no_check_fields=[]):
    # Suppress warnings to work with BAM files without index file. 
    # https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    # save = pysam.set_verbosity(0)
    # bam_file = pysam.AlignmentFile(bam_path, "rb")
    # pysam.set_verbosity(save)

    fields_to_check = [field for field in list(map(int, Fields)) if field not in no_check_fields]

    gbam_file = get_reader(gbam_path, get_parsing_tmpl(fields_to_check)) 
    from gbam_tools import GbamRecord

    while True:
        cur_gbam = gbam_file.next_rec()
        print(cur_gbam)
        # cur_bam = next(bam_file, None)

        # if cur_gbam == None or cur_bam == None:
        #     # Assert there is no records left
        #     assert(cur_gbam == cur_bam)
        #     break
        
        # for field in fields_to_check:
        #     if field == Fields.RAWSEQUENCE:
        #         assert(cur_bam.query_sequence == cur_gbam.seq)
        #     if field == Fields.REFID:
        #         assert(cur_bam.reference_id == cur_gbam.refid)
        #     if field == Fields.MAPQ:
        #         assert(cur_bam.mapping_quality == cur_gbam.mapq)
        #     if field == Fields.POS:
        #         assert(cur_bam.reference_start == cur_gbam.pos)
        #     if field == Fields.RAWQUAL:
        #         assert(cur_bam.query_qualities == cur_gbam.qual)
        #     if field == Fields.RAWCIGAR:
        #         assert(cur_bam.cigarstring == cur_gbam.cigar)
        #     if field == Fields.READNAME:
        #         assert(cur_bam.reference_name == cur_gbam.read_name)
            
        # else:
        #     assert(cur_gbam.mapq == None)
        # if Fields.POS not in no_check_fields:
            
        # else:
        #     assert(cur_gbam.pos == None)
                
def is_valid_file(parser, arg):
    import os.path

    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Tool for creating and reading GBAM files",
    )

    subparsers = parser.add_subparsers(
        help='Choose one of the following options.', dest="command")
    subparsers.required = True

    parser_convert = subparsers.add_parser(
        'conv', help='Convert BAM file to GBAM file.')
    parser_convert.add_argument("-i", "--input_path", help="BAM file to convert.", dest="input_path",
                                type=lambda x: is_valid_file(parser_convert, x), required=True)
    parser_convert.add_argument(
        "-o", "--output_file", type=str, help="GBAM file name to write out.", required=True)

    # parser_read = subparsers.add_parser('read', help='Read GBAM file.')
    # parser_read.add_argument("-i", "--input_path", help="GBAM file to read.", dest="input_path",
    #                          type=lambda x: is_valid_file(parser_convert, x), required=True)
    # parser_read.add_argument('--no_mapq', dest='mapq', action='store_true',
    #                          help='Avoid reading MAPQ field.')
    # parser_read.add_argument('--no_pos', dest='pos', action='store_true',
    #                          help='Avoid reading POS field.')

    parser_test = subparsers.add_parser('test', help='Test GBAM tools.')
    parser_test.add_argument("-i", "--input_path", help="GBAM file to perform tests on.", dest="input_path",
                             type=lambda x: is_valid_file(parser_convert, x), required=True)

    import sys

    args = parser.parse_args(sys.argv[1:])

    if args.command == 'conv':
        convert(args.input_path, args.output_file)
    # elif args.command == 'read':
    #     args = parser_read.parse_args(sys.argv[1:])
    #     read(args)
    elif args.command == 'test':
        test(args)
