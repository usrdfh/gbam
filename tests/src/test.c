#include "../../c_api/inc/api.h"
#include <dlfcn.h>
#include <assert.h>
#include <unistd.h>

int main(int argc, char **argv){
    void* reader = get_gbam_reader(argv[1], DISABLE_TAGS | DISABLE_SEQ);
    sam_hdr_t* hdr = get_header(reader);

    int amount = 0;
    bam1_t* b = get_empty_bam1_t();
    while(1){
        if(fill_next_record(reader, b) == -1){
            break;
        }
        // kstring_t s = {0,0,NULL};
        // int res = sam_format1(hdr, rec, &s);
        // printf("%s\ n", s.s);
        // free(s.s);
        amount++;
    }
    // free_bam_record(b);
    // for(int i = 0; i < get_records_num(reader); i++){
       
    //     bam1_t* rec = get_bam_record_fast(reader, i);
  
    //     amount++;
        
       

    // }
    printf("%d", amount);
    free_header(hdr);
}