SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
gcc -o $SCRIPT_DIR/../target/release/ffi_test.o $SCRIPT_DIR/src/test.c -L$SCRIPT_DIR/../target/release/ -lgbam_tools_cffi -lhts