#ifndef MMLIB_H
#define MMLIB_H

#include <stddef.h>
#include <inttypes.h>

#define HEADMASK 0x80000000

typedef struct {
	size_t n;
	size_t m;
	uint32_t *a;
} uint32_v_t;

/*
seed has been added following 
https://github.com/medvedevgroup/minimizer-jaccard-estimator/blob/main/minimap2_hash_uncompiled.py
version in order to provide independent trials, if needed
*/
static inline uint64_t hash64(uint64_t seed, uint64_t key, uint64_t mask)
{
	key = (key + seed) & mask;
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

int mm_get_pos_pool(void *km, const char *seq, size_t slen, uint8_t m, uint16_t w, uint64_t seed, uint32_v_t *mm_pos);

int mm_get_pos(const char *seq, size_t slen, uint8_t m, uint16_t w, uint64_t seed, uint32_v_t *mm_pos);

int sync_get_pos_pool(void *km, const char *seq, size_t slen, uint8_t k, uint8_t s, uint64_t seed, uint32_v_t *sync_pos);

int sync_get_pos(const char *seq, size_t slen, uint8_t k, uint8_t s, uint64_t seed, uint32_v_t *sync_pos);

#endif/*MMLIB_H*/
