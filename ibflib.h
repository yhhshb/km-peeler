#ifndef IBFLIB_H
#define IBFLIB_H

#include "constants.h"

#define MAXCELLLEN (WSIZE * 8)

/* FIXME [possible bug] Useless because lengths must use + instead of ^ */
#if (MAXCELLLEN <= 0xFF)
#define keysum_len_t uint8_t
#define hton_len(v) (hton8(v))
#define ntoh_len(v) (ntoh8(v))
#define format_len "hhu"
#elif (MAXCELLLEN <= 0xFFFF)
#define keysum_len_t uint16_t
#define hton_len(v) (hton16(v))
#define ntoh_len(v) (ntoh16(v))
#define format_len "hu"
#elif (MAXCELLLEN <= 0xFFFFFFFF)
#define keysum_len_t uint32_t
#define hton_len(v) (hton32(v))
#define ntoh_len(v) (ntoh32(v))
#define format_len "u"
#elif (MAXCELLLEN <= 0xFFFFFFFFFFFFFFFF)
#define keysum_len_t uint64_t
#define hton_len(v) (hton64(v))
#define ntoh_len(v) (ntoh64(v))
#define format_len "llu"
#else
#define hton_len(v) (unsupported_length_type(v))
#define ntoh_len(v) (unsupported_length_type(v))
#endif

#if defined(RPOS)
/*do nothing*/
#elif (MAXSEQLEN <= 0xFF)
#define position_t uint8_t
#define hton_pos(v) (hton8(v))
#define ntoh_pos(v) (ntoh8(v))
#define format_pos "hhu"
#elif (MAXSEQLEN <= 0xFFFF)
#define position_t uint16_t
#define hton_pos(v) (hton16(v))
#define ntoh_pos(v) (ntoh16(v))
#define format_pos "hu"
#elif (MAXSEQLEN <= 0xFFFFFFFF)
#define position_t uint32_t
#define hton_pos(v) (hton32(v))
#define ntoh_pos(v) (ntoh32(v))
#define format_pos "u"
#elif (MAXSEQLEN <= 0xFFFFFFFFFFFFFFFF)
#define position_t uint64_t
#define hton_pos(v) (hton64(v))
#define ntoh_pos(v) (ntoh64(v))
#define format_pos "llu"
#else
#define hton_pos(v) (unsupported_length_type(v))
#define hton_pos(v) (unsupported_length_type(v))
#endif

typedef struct {
    uint64_t ls64b;
    uint64_t ms64b;
} hash_t;

typedef struct {
    uint32_t seed;/*seed for a block (constant)*/
    hash_t hash;/*result of the hashing function (changes when a sequence is hashed)*/
} hash_gen_t;

typedef struct {
    int64_t counter;/*<POSSIBLE SOURCE OF ERRORS: changed from unsigned to signed because of symmetry. If bugs use a defined threshold = 2^63*/
    uint8_t keysum[WSIZE];
    #ifndef GLEN/*if all keys are the same length, then store it inside the ibf itself, saving space in the buckets*/
    keysum_len_t key_len;
    #endif
    #ifndef RPOS
    position_t position;
    #endif
} bucket_t;

typedef struct {
    /*uint32_t seed;*/
    uint8_t repetitions;/* number of hashes/blocks */
    float epsilon;/*approximation factor*/
    uint64_t chunk_size;/*depends on r, and the expected number of differences (+ the approx factor to augment the prob. of success)*/
    #ifdef GLEN/*if all keys are the same length, it is stored here and not into each bucket*/
    keysum_len_t key_len;
    #endif
    bucket_t *data;/*the sketch itself*/
    hash_gen_t *seres;/* seeds + results for each block */
} ibf_t;

int ibf_sketch_init(unsigned int root_seed, unsigned char r, float epsilon, unsigned int n, ibf_t *const sketch);

int ibf_sketch_copy(ibf_t const *const source, ibf_t *const dest);

int ibf_sketch_destroy(ibf_t *const sketch);

int ibf_sketch_store(char const * const path, ibf_t const * const sketch);

int ibf_sketch_load(char const * const path, ibf_t * const sketch);

int ibf_sketch_print(ibf_t const *const sketch, FILE *const strm);

int ibf_sketch_dump(ibf_t const *const sketch, FILE *const strm);

int ibf_sketch_diff(ibf_t const *const a, ibf_t const *const b, ibf_t *const result);

int ibf_insert_seq(void const *const seq, int start, int end, ibf_t *const sketch, uint8_t *const buffer, uint64_t buffer_len);

int ibf_delete_seq(void const *const seq, int start, int end, ibf_t *const sketch, uint8_t *const buffer, uint64_t buffer_len);

int ibf_list_seq(ibf_t *const sketch, void (*output_bucket)(bucket_t const *const, char, void*), void *iostruct);/*DESTRUCTIVE OPERATION, make copy of sketch if needed*/

int ibf_count_seq(ibf_t const *const sketch, unsigned long *const count);

int ibf_buffer_init(uint8_t** const buffer);

int ibf_buffer_destroy(uint8_t** const buffer);

#endif/*IBFLIB*/
