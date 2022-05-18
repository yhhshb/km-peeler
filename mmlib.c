#include <string.h>
#include "kvec2.h"
#include "err.h"
#include "mmlib.h"
#include "compile_options.h"
#include "constants.h"

/*#include <stdio.h>*/
#include <assert.h>

typedef struct {
    uint64_t hash;
    uint32_t pos;
} mm_t;

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param seq    DNA sequence
 * @param slen   sequence length
 * @param k      k-mer size
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param mm_pos minimizer positions in the string
 *               where a position is the index of the first base of a minimizer,
 *               Callers may want to set "mm_pos->n = 0"; otherwise results are appended to mm_pos
 */
int mm_get_pos_pool(void *km, const char *seq, size_t slen, uint8_t m, uint16_t w, uint64_t seed, uint32_v_t *mm_pos)
{
	int c;
	unsigned char z;
	uint64_t shift, mask, kmer[2];
	mm_t buf[256], min, info;
	size_t i, j, pns_last_inv_base, buf_pos, min_pos;

	assert(seq != NULL);
	assert(mm_pos != NULL);
	/*
	assert(slen <= MAXSEQLEN);
	assert(m <= MMML);
	assert(w <= MWW);
	*/

	kv_resize(uint32_t, km, *mm_pos, mm_pos->n + (slen / w));

	shift = 2 * (m - 1);
	mask = (1ULL<<2 * m) - 1;
	memset(kmer, 0, sizeof kmer);
	memset(buf, 0xFF, w * 16);
	min.hash = UINT64_MAX, min.pos = UINT32_MAX;
	pns_last_inv_base = buf_pos = min_pos = 0;
	z = 0;

	for (i = 0; i < slen; ++i) {
		c = seq_nt4_table[(uint8_t)seq[i]];
		//fprintf(stderr, "reading base: %d, (seq[%lu] = %c\n", c, i, seq[i]);
		info.hash = UINT64_MAX, info.pos = UINT32_MAX;
		if (c < 4) {
			kmer[0] = (kmer[0] << 2 | c) & mask;          /* forward k-mer */
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift; /* reverse k-mer */
			if (kmer[0] != kmer[1]) z = kmer[0] < kmer[1] ? 0 : 1; /* strand, if symmetric k-mer then use previous strand */
			++pns_last_inv_base;
			if (pns_last_inv_base >= m) {
				info.hash = hash64(seed, kmer[z], mask);
				info.pos = (uint32_t)(i-m+1);/* +1 because i starts from 0 but m is the actual length*/
			}
		} else {
			pns_last_inv_base = 0;
		}
		//fprintf(stderr, "buf[%lu] = (%llu, %u)\n", buf_pos, buf[buf_pos].hash, buf[buf_pos].pos);
		if (pns_last_inv_base != 0) {
			//fprintf(stderr, "info = (%llu, %u) -> buf[%lu]\n", info.hash, info.pos, buf_pos);
			buf[buf_pos] = info; /* need to do this here as appropriate buf_pos and buf[buf_pos] are needed below */
		}
		if (pns_last_inv_base == w + m - 1 && min.hash != UINT64_MAX) { /* special case for the first window - because identical k-mers are not stored yet */
			//fprintf(stderr, "first window pns_last_inv_base = %lu\n", pns_last_inv_base);
			for (j = buf_pos + 1; j < w; ++j)
				if (min.hash == buf[j].hash && buf[j].pos != min.pos) kv_push(uint32_t, km, *mm_pos, buf[j].pos);
			for (j = 0; j < buf_pos; ++j)
				if (min.hash == buf[j].hash && buf[j].pos != min.pos) kv_push(uint32_t, km, *mm_pos, buf[j].pos);
		}
		if (info.hash <= min.hash) { /* a new minimum; then write the old min */
			//fprintf(stderr, "info.hash <= min.hash ((%llu, %u), (%llu, %u))", info.hash, info.pos, min.hash, min.pos);
			if (pns_last_inv_base >= w + m && min.hash != UINT64_MAX) {
				//fprintf(stderr, ", writing min.pos = %u", min.pos);
				kv_push(uint32_t, km, *mm_pos, min.pos);
			}
			//fprintf(stderr, "\n");
			min = info;
			min_pos = buf_pos;
		} else if (buf_pos == min_pos) { /* old min has moved outside the window */
			//fprintf(stderr, "buf_pos == min_pos (%lu == %lu)\n", buf_pos, min_pos);
			if ((pns_last_inv_base ==0 || pns_last_inv_base >= w + m - 1) && min.hash != UINT64_MAX) kv_push(uint32_t, km, *mm_pos, min.pos);
			for (j = buf_pos + 1, min.hash = UINT64_MAX; j < w; ++j) {/* the two loops are necessary when there are identical k-mers */
				if (min.hash >= buf[j].hash) {
					min = buf[j];
					min_pos = j;/* >= is important s.t. min is always the closest k-mer */
				}
			}
			for (j = 0; j <= buf_pos; ++j) {
				if (min.hash >= buf[j].hash && buf[j].hash != UINT64_MAX) {
					min = buf[j];
					min_pos = j;
				}
			}
			if (pns_last_inv_base >= w + m - 1 && min.hash != UINT64_MAX) {/* write identical k-mers */
				/* these two loops make sure the output is sorted */
				for (j = buf_pos + 1; j < w; ++j) 
					if (min.hash == buf[j].hash && min.pos != buf[j].pos) kv_push(uint32_t, km, *mm_pos, buf[j].pos);
				for (j = 0; j <= buf_pos; ++j) 
					if (min.hash == buf[j].hash && min.pos != buf[j].pos) kv_push(uint32_t, km, *mm_pos, buf[j].pos);
			}
		}
		if (w == 0 || ++buf_pos == w) buf_pos = 0;
		//fprintf(stderr, "buf_pos = %lu\n", buf_pos);
	}
	if (min.hash != UINT64_MAX)
		kv_push(uint32_t, km, *mm_pos, min.pos);
    return NO_ERROR;
}

int mm_get_pos(const char *seq, size_t slen, uint8_t m, uint16_t w, uint64_t seed, uint32_v_t *mm_pos) {
	return mm_get_pos_pool(NULL, seq, slen, m, w, seed, mm_pos);
}

int sync_get_pos_pool(void *km, const char *seq, size_t slen, uint8_t k, uint8_t s, uint64_t seed, uint32_v_t *sync_pos) {
	int c;
	unsigned char z;
	uint64_t shift, mask, smer[2];
	mm_t buf[256], info;
	size_t i, j, buf_pos, min_pos, pns_last_inv_base;
	const size_t w = k-s+1;

	assert(seq != NULL);
	assert(sync_pos != NULL);
	
	shift = 2 * (s - 1);
	mask = (1ULL<<2 * s) - 1;
	memset(smer, 0, sizeof smer);
	memset(buf, 0xFF, k * 16);
	//min.hash = UINT64_MAX, min.pos = UINT32_MAX;
	buf_pos = pns_last_inv_base = 0;
	min_pos = UINT64_MAX;
	z = 0;

	for(i = 0; i < slen; ++i) {
		c = seq_nt4_table[(uint8_t)seq[i]];
		/*fprintf(stderr, "reading base: %c (seq[%lu])\n", seq[i], i);*/
		info.hash = UINT64_MAX, info.pos = UINT32_MAX;
		if (c < 4) {
			smer[0] = (smer[0] << 2 | c) & mask;          /* forward s-mer */
			smer[1] = (smer[1] >> 2) | (3ULL^c) << shift; /* reverse s-mer */
			if (smer[0] != smer[1]) z = smer[0] < smer[1] ? 0 : 1; /* strand, it does not matter and hash(smer[0]) == hash(smer[1]) if smer[0] == smer[1] */
			++pns_last_inv_base;
			if (pns_last_inv_base >= s) {
				info.hash = hash64(seed, smer[z], mask);
				info.pos = (uint32_t)(i-s+1);
			}
		} else {
			pns_last_inv_base = 0;
		}
		/*fprintf(stderr, "pns = %lu | min_pos = %lu | buf_pos = %lu\n", pns_last_inv_base, min_pos, buf_pos);*/
		buf[buf_pos] = info;
		if (buf_pos == min_pos) min_pos = UINT64_MAX;/*old minimum overwritten by new s-mer*/
		if (pns_last_inv_base >= k) {
			if (min_pos == UINT64_MAX) {/*find new minimum if old out of window*/
				min_pos = buf_pos;
				for(j = buf_pos + 1; j < k - s + 1; ++j) if (buf[j].hash < buf[min_pos].hash) min_pos = j;
				for(j = 0; j <= buf_pos; ++j) if (buf[j].hash < buf[min_pos].hash) min_pos = j;
				/*ATTENTION HERE          |  -- Before it was < (strictly inferior).
				                          |  	which gives different results depending on if
										  |		k-mer set (basically a fasta with all sequences of length k) or 
										  ▼		the original string were given as input.*/
			} else if (buf[buf_pos].hash <= buf[min_pos].hash) {/*update minimum if still inside window*/
				/*fprintf(stderr, "min_pos = buf_pos = %lu\n", buf_pos);*/
				min_pos = buf_pos;
			}
			if (min_pos == ((buf_pos + 1) % w)) {/*syncmers with min at the beginning*/
				/*fprintf(stderr, "syncmer found (min_pos = %lu | buf_pos = %lu)\n", min_pos, buf_pos);*/
				kv_push(uint32_t, km, *sync_pos, buf[min_pos].pos);
				/*
				fprintf(stderr, "minimum at the beginning:\n");
				for(j = buf_pos + 1; j < k - s + 1; ++j) fprintf(stderr, "%llu,", buf[j].hash);
				for(j = 0; j <= buf_pos; ++j) fprintf(stderr, "%llu,", buf[j].hash);
				fprintf(stderr, "\n");
				*/
			} else if (min_pos == buf_pos) {/*syncmer with min at the end*/
				kv_push(uint32_t, km, *sync_pos, buf[min_pos].pos - k + s);
				/*
				fprintf(stderr, "minimum at the end:\n");
				for(j = buf_pos + 1; j < k - s + 1; ++j) fprintf(stderr, "%llu,", buf[j].hash);
				for(j = 0; j <= buf_pos; ++j) fprintf(stderr, "%llu,", buf[j].hash);
				fprintf(stderr, "\n");
				*/
			}
		}
		if (++buf_pos == w || k == 0) buf_pos = 0;
	}
	return NO_ERROR;
}

int sync_get_pos(const char *seq, size_t slen, uint8_t k, uint8_t s, uint64_t seed, uint32_v_t *sync_pos) {
	return sync_get_pos_pool(NULL, seq, slen, k, s, seed, sync_pos);
}
