#ifndef MINHASH_H
#define MINHASH_H

#include <stdint.h>

typedef struct {
    uint32_t seed;
    uint64_t size;
    uint8_t hash_width;
    uint8_t *hashes;
} minhash_t;

int minhash_sketch_init(uint32_t seed, uint64_t s, uint8_t hwidth, minhash_t *const sketch);

int minhash_sketch_to_stream(FILE* ostrm, minhash_t const * const sketch);

int minhash_sketch_store(char const *const path, minhash_t const *const sketch);

int minhash_sketch_from_stream(FILE* istrm, minhash_t *const sketch);

int minhash_sketch_load(char const *const path, minhash_t *const sketch);

int minhash_print(FILE* out, minhash_t const * const sketch);

int minhash_sketch_insert(char const *const kmer, uint8_t k, uint64_t s, minhash_t *const sketch);

int minhash_compare(minhash_t const *const s1, minhash_t const *const s2, double *jaccard);

int minhash_sketch_destroy(minhash_t *const sketch);

#endif/*MINHASH_H*/
