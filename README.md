# GBAM

This package includes tools for creating and manipulating GBAM files. GBAM is a new column oriented file format for storing binary alignment data. The way it stores information allows to skip unnecessary operations which are inevitable when using regular BAM files.

# Installation

## Install binary

```shell
git clone https://github.com/NickRoz1/gbam
cd gbam
cargo build --release
./target/release/gbam_binary --help
```

You may need to install following packages:
```shell
sudo apt-get install liblzma-dev
sudo apt-get install libbz2-dev
sudo apt-get install libhts-dev
```

# Usage

### Examples
```shell
# Simply convert
time ./target/release/gbam_binary -c test.bam -o test.gbam

# Sort before writing (sort by reference and coordinates (other sort predicates are available, but not implemented in CLI currently))
time ./target/release/gbam_binary -c -s 1gb.bam -o 1gb.sorted.gbam --sort-temp-mode [lz4_file|file|lz4_ram|ram]

# Collect flag statistics
time ./target/release/gbam_binary --flagstat test.gbam

# Calculate read depth (only on sorted files)
time ./target/release/gbam_binary --depth test.sorted.gbam --thread-num 4 > depth_test.txt

# Calculate read depth (only on sorted files) and create bed regions depth gzip file
time ./target/release/gbam_binary --depth test.sorted.gbam --thread-num 4 -o test_data/depth_test.bed.gz
```

### To run pytests
```shell
# Run all tests
pytest

# To run all correctness test on a specific file path
pytest --use-custom-file=/home/test/testing_big_file/test.bam

# Test sort on user provided file (slow)
pytest -k test_sort --bam-file-path=/home/test/testing_big_file/test.bam

# Test conversion bam to gbam to bam (preserving correct data) on user provided file
pytest -k test_bam_to_gbam_to_bam --bam-file-path=/home/test/testing_big_file/test.bam
```

### Benchmarks

See [benchmarks](benchmarking/benchmark_result.md)

Data is available for [download](http://hypervolu.me/~guarracino/gbam_injection/).

# Development

## GNU Guix - this is not tested.

GNU Guix provides a full environment for development.  See
[.guix-build](./.guix-build)
