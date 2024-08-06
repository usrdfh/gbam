// api.h
#ifndef API_H
#define API_H

#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

static const uint64_t DISABLE_READNAME = 1<<0;
static const uint64_t DISABLE_CIGAR = 1<<1;
static const uint64_t DISABLE_SEQ = 1<<2;
static const uint64_t DISABLE_QUAL = 1<<3;
static const uint64_t DISABLE_TAGS = 1<<4;

void* get_gbam_reader(const char* gbam_path, uint64_u fields_to_not_set);
bam1_t get_bam_record(void* reader, uint64_u rec_num);
bam1_t* get_empty_bam1_t();
bam1_t* get_bam_record_fast(void* reader, uint64_u rec_num);
int32_t fill_next_record(void*  reader, bam1_t* rec);

sam_hdr_t* get_header(void* reader);
uint64_t get_records_num(void* reader);

void free_bam_record(bam1_t* rec);
void free_gbam_reader(void* reader);
void free_header(sam_hdr_t* hdr);

#ifdef __cplusplus
}
#endif

#endif // API_H