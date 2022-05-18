#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "err.h"
#include "constants.h"
#include "endian_fixer.h"
#include "murmur3.h"
#include "minHash.h"

#include <assert.h>

int minhash_sketch_init(uint32_t seed, uint64_t s, uint8_t hwidth, minhash_t *const sketch) {
    uint64_t byte_size;
    assert(sketch);
    sketch->seed = seed;
    sketch->size = 0;
    sketch->hash_width = CEILING(hwidth, 8);
    byte_size = (s + 1) * sketch->hash_width;
    if((sketch->hashes = (uint8_t*)malloc(byte_size)) == NULL) return ERR_ALLOC; /* +1 is to avoid overflows with memmove */
    return NO_ERROR;
}

int minhash_sketch_to_stream(FILE* ostrm, minhash_t const *const sketch) {
    uint32_t buffer32;
    uint64_t buffer64;
    assert(ostrm);
    assert(sketch);
    buffer32 = hton32(sketch->seed);
    if (fwrite(&buffer32, sizeof buffer32, 1, ostrm) != 1) return ERR_IO;
    buffer64 = hton64(sketch->size);
    if (fwrite(&buffer64, sizeof buffer64, 1, ostrm) != 1) return ERR_IO;
    if (fwrite(&sketch->hash_width, sizeof sketch->hash_width, 1, ostrm) != 1) return ERR_IO;
    if (fwrite(sketch->hashes, sketch->hash_width, sketch->size, ostrm) != sketch->size) return ERR_IO;
    return NO_ERROR;
}

int minhash_sketch_store(char const *const path, minhash_t const *const sketch) {
    FILE* out;
    int err;
    assert(path);
    assert(sketch);
    if((out = fopen(path, "wb")) == NULL) return ERR_IO;
    err = minhash_sketch_to_stream(out, sketch);
    fclose(out);
    return err;
}

int minhash_sketch_from_stream(FILE* istrm, minhash_t *const sketch) {
    assert(istrm);
    assert(sketch);
    if (fread(&sketch->seed, sizeof sketch->seed, 1, istrm) != 1) return ERR_IO;
    if (fread(&sketch->size, sizeof sketch->size, 1, istrm) != 1) return ERR_IO;
    if (fread(&sketch->hash_width, sizeof sketch->hash_width, 1, istrm) != 1) return ERR_IO;
    sketch->seed = ntoh32(sketch->seed);
    sketch->size = hton64(sketch->size);
    if ((sketch->hashes = malloc(sketch->hash_width * sketch->size)) == NULL) return ERR_ALLOC;
    if (fread(sketch->hashes, sketch->hash_width, sketch->size, istrm) != sketch->size) return ERR_IO;
    return NO_ERROR;
}

int minhash_sketch_load(char const *const path, minhash_t *const sketch) {
    FILE* in;
    int err;
    assert(path);
    assert(sketch);
    if ((in = fopen(path, "rb")) == NULL) return ERR_IO;
    err = minhash_sketch_from_stream(in, sketch);
    fclose(in);
    return err;
}

static inline void* at(minhash_t const *const sketch, uint64_t idx) {
    return &sketch->hashes[idx * sketch->hash_width];
}

static inline void print_hex(FILE* out, uint8_t const *const p, uint64_t length) {
    uint64_t i;
    for(i = 0; i < length; ++i) fprintf(out, "%X ", p[i]);
}

static inline void print_dec(FILE* out, uint8_t const * const p) {
    fprintf(out, "%llu", *((uint64_t*)p));
}

int minhash_print(FILE* out, minhash_t const * const sketch) {
    uint64_t i;
    assert(out);
    assert(sketch);
    for(i = 0; i < sketch->size; ++i) {
        if(sketch->hash_width == 8) print_dec(out, at(sketch, i));
        else print_hex(out, at(sketch, i), sketch->hash_width);
        fprintf(out, " : %llu\n", i);
    }
    return NO_ERROR;
}

int minhash_sketch_insert(char const *const kmer, uint8_t k, uint64_t s, minhash_t *const sketch) {
    uint64_t i;
    int cmp;
    uint8_t buffer[16];
    assert(kmer);
    assert(sketch);
    MurmurHash3_x64_128(kmer, k, sketch->seed, buffer);
    // fprintf(stderr, "%.*s = ", (int)k, kmer);
    // print_hex(buffer, 16);
    // fprintf(stderr, "\n");
    cmp = -1;
    for(i = sketch->size - 1; i != UINT64_MAX && (cmp = memcmp(buffer, at(sketch, i), sketch->hash_width)) < 0; --i) {}/* find insertion index */
    //fprintf(stderr, "sketch size = %llu, cmp = %d\n", sketch->size, cmp);
    if (cmp == 0) {
        // fprintf(stderr, "collision kmer = %.*s -> %llu\n", k, kmer, ((uint64_t*)buffer)[0]);
        return NO_ERROR;
    }
    if (++i < s) {/* shift and insert new hash if new minimum or if sketch not yet full */
        // fprintf(stderr, "[Insert] ");
        // print_hex(stderr, buffer, sketch->hash_width);
        // fprintf(stderr, " at sketch[%llu] =  ", i);
        // print_hex(stderr, at(sketch, i), sketch->hash_width);
        // fprintf(stderr, "\n");
        memmove(at(sketch, i+1), at(sketch, i), (s - i) * sketch->hash_width);
        memcpy(at(sketch, i), buffer, sketch->hash_width);
        if (sketch->size < s) ++sketch->size;
    }
    return NO_ERROR;
}

int minhash_compare(minhash_t const *const s1, minhash_t const *const s2, double *jaccard) {
    uint64_t num, den, s;
    uint64_t i, j;
    assert(s1);
    assert(s2);
    assert(jaccard);
    if (s1->hash_width != s2->hash_width || s1->seed != s2->seed) return ERR_INCOMPATIBLE;
    if (s1->size < s2->size) s = s1->size;
    else s = s2->size;
    num = den = 0;
    i = j = 0;
    while(i < s1->size && j < s2->size && den < s) {
        int cmp = memcmp(at(s1, i), at(s2, j), s1->hash_width);
        if (cmp < 0) ++i;
        else if (cmp > 0) ++j;
        else {++i; ++j; ++num;}
        ++den;
    }
    *jaccard = ((double)num) / (double)den;
    fprintf(stderr, "%f = %llu/%llu\n", *jaccard, num, den);
    return NO_ERROR;
}

int minhash_sketch_destroy(minhash_t *const sketch) {
    assert(sketch);
    if (sketch->hashes) {
        free(sketch->hashes);
        sketch->hashes = NULL;
    }
    return NO_ERROR;
}
