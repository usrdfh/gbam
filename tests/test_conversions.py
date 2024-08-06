from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import pytest 
import gzip
import shutil
import os
import io

with_depth = pytest.mark.skipif("not config.getoption('with_depth')")

cur_file_path = Path(__file__).parent.absolute()

test_data_folder = cur_file_path.parent/"test_data"
binary_path = cur_file_path.parent/"target"/"release"/"gbam_binary"

bam_file_path = test_data_folder/"little.bam"

gbam_file = None
gbam_file_sorted = None
bam_file_sorted_path = None

@pytest.fixture(scope='module', autouse=True)
def gen_gbam_file(request):
    global gbam_file, gbam_file_sorted, bam_file_sorted_path, bam_file_path
    if request.config.getoption("--use-custom-file") is not None:
        bam_file_path = Path(request.config.getoption("--use-custom-file"))
    gbam_file = NamedTemporaryFile()
    # Simply convert to GBAM.
    subprocess.run([binary_path, bam_file_path, "-c", "-o", gbam_file.name]) 
    gbam_file_sorted = NamedTemporaryFile()
    # Sort and convert to GBAM.
    subprocess.run([binary_path, bam_file_path, "-c", "-s", "-o", gbam_file_sorted.name, "--index-sort"]) 
    bam_file_sorted_path = NamedTemporaryFile(suffix=".bam")
    subprocess.run(["samtools", "sort", "-@ 8", bam_file_path, "-o", bam_file_sorted_path.name]) 

def compare_bam_files(original, result):
    out_of_original_view = NamedTemporaryFile()
    out_of_result_view = NamedTemporaryFile()
    
    subprocess.check_call(["samtools", "view", "-@ 8", original, '-o', out_of_original_view.name], stderr=subprocess.STDOUT)
    subprocess.check_call(["samtools", "view", "-@ 8", result, '-o', out_of_result_view.name], stderr=subprocess.STDOUT)

    byte_file_comparison(out_of_original_view.name, out_of_result_view.name)
   

def byte_file_comparison(original_path, result_path):
    assert(os.path.getsize(result_path) > 0)
    assert(os.path.getsize(result_path) == os.path.getsize(original_path))
    
    with io.open(original_path, 'rb') as from_samtools, io.open(result_path, 'rb') as from_gbam:
        while True:
            bam_byte = from_samtools.read(4096)
            gbam_byte = from_gbam.read(4096)

            if not bam_byte:
                break

            assert(bam_byte == gbam_byte)

def test_bam_to_gbam_to_bam(request):
    bam_file_from_gbam = NamedTemporaryFile(suffix=".bam")
    test_bam_file_path = bam_file_path
    cli_bam_path = request.config.getoption("--bam-file-path")
    
    if cli_bam_path is None:
        subprocess.run([binary_path, "--convert-to-bam", gbam_file.name, "-o", bam_file_from_gbam.name]) 
    else:
        cli_bam_path = Path(cli_bam_path)
        temp_gbam = NamedTemporaryFile()
        subprocess.run([binary_path, cli_bam_path.as_posix(), "-c", "-o", temp_gbam.name]) 
        subprocess.run([binary_path, "--convert-to-bam", temp_gbam.name, "-o", bam_file_from_gbam.name]) 
        test_bam_file_path = cli_bam_path
        del temp_gbam

    compare_bam_files(test_bam_file_path.as_posix(), bam_file_from_gbam.name)

def test_flagstat():
    view_of_original = subprocess.check_output(["samtools", "flagstat", str(bam_file_path)], stderr=subprocess.STDOUT)
    view_of_result = subprocess.check_output([binary_path, "--flagstat", gbam_file.name], stderr=subprocess.STDOUT)

    assert(len(view_of_original) > 0)
    assert(view_of_original == view_of_result)

def test_sort(request):
    gbam_sorted_results = NamedTemporaryFile(suffix=".bam")
    samtools_sorted_results = bam_file_sorted_path
    cli_bam_path = request.config.getoption("--bam-file-path")
    
    if cli_bam_path is None:
        gbam_res, samtools_res = generate_views_for_gbam_and_bam_files(gbam_file_sorted, gbam_file_sorted.name+".gbai", bam_file_sorted_path.name)
        byte_file_comparison(gbam_res.name, samtools_res.name)
        return
    
    cli_bam_path = Path(cli_bam_path)
    temporary_gbam = NamedTemporaryFile()
    subprocess.run([binary_path, cli_bam_path.as_posix(), "-c", "-s", "-o", temporary_gbam.name, "--sort-temp-mode", "lz4_ram"]) 
    subprocess.run([binary_path, "--convert-to-bam", temporary_gbam.name, "-o", gbam_sorted_results.name]) 
    samtools_sorted_results = NamedTemporaryFile(suffix=".bam")
    subprocess.run(["samtools", "sort", cli_bam_path.as_posix(), "-o", samtools_sorted_results.name]) 

    compare_bam_files(samtools_sorted_results.name, gbam_sorted_results.name)
    
# Testing against mosdepth.
@with_depth
def test_depth():    
    temp_dir = TemporaryDirectory()
    dir_path = Path(temp_dir.name)

    # Mosdepth won't work without index and won't accept absolute paths. Copy BAM file and index file into mosdepth temp directory.
    shutil.copy(bam_file_sorted_path.name, dir_path)
    bam_file_name = bam_file_sorted_path.name.split("/")[-1]
    bam_file_path = dir_path/bam_file_name
    index_file_path = dir_path/(bam_file_name+".bai")
    subprocess.check_call(['samtools', 'index', bam_file_sorted_path.name, '-o', index_file_path.as_posix()])
    
    mosdepth_prefix = "testing_depth_pytest"
    mosdepth_suffix = ".per-base.bed.gz"
    mosdepth_file = (dir_path/(mosdepth_prefix+mosdepth_suffix)).as_posix()

    subprocess.check_call(['mosdepth', '-x', mosdepth_prefix, bam_file_path.name], cwd=temp_dir.name)
    bed_gz = NamedTemporaryFile()
    subprocess.check_call([binary_path, gbam_file_sorted.name, '-d', '-o', bed_gz.name, '--index-file', gbam_file_sorted.name + '.gbai'])
    
    with gzip.open(bed_gz, "rb") as gbam_res, gzip.open(mosdepth_file, "rb") as mosdepth_res:
        # Read the decompressed data from both files
        decompressed_gbam_res = gbam_res.read()
        decompressed_mosdepth_res= mosdepth_res.read()

    assert(len(decompressed_gbam_res) > 0)
    assert(decompressed_gbam_res == decompressed_mosdepth_res)

def generate_views_for_gbam_and_bam_files(gbam_input, gbam_index, bam_input):
    gbam_results = NamedTemporaryFile()
    samtools_results = NamedTemporaryFile()

    index_option = ""
    if gbam_index is not None:
        index_option = f"--index-file {gbam_index}"

    subprocess.check_output([f"{binary_path} -v {gbam_input.name} {index_option} | samtools view > {gbam_results.name}"], shell=True) 
    subprocess.check_output([f"samtools view {bam_input} > {samtools_results.name}"], shell=True) 

    return (gbam_results, samtools_results)

def test_view(request):
    gbam_results, samtools_results = generate_views_for_gbam_and_bam_files(gbam_file, None, bam_file_path)
    byte_file_comparison(samtools_results.name, gbam_results.name)

def test_cffi():
    gbam_results = NamedTemporaryFile()
    samtools_results = NamedTemporaryFile()

    # Build the .c test
    subprocess.check_output([f"bash {cur_file_path.parent}/tests/build_test.sh"], shell=True) 

    subprocess.check_output([f"samtools view {bam_file_path} > {samtools_results.name}"], shell=True) 

    test_bin = cur_file_path.parent/"target"/"release"/"ffi_test.o"
    subprocess.check_output([f"{test_bin} {gbam_file.name} > {gbam_results.name}"], shell=True) 

    byte_file_comparison(samtools_results.name, gbam_results.name)
