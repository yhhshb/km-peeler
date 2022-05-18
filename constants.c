#include <stddef.h>
#include "constants.h"
#include "err.h"

#include <assert.h>

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char seq_nt4_inv_table[5] = {'A','C','G','T','N'};

int pack2bit(const char *seq, unsigned char len, unsigned char *out) {
	int i;
	unsigned char c;
	assert(seq != NULL);
	assert(out != NULL);
	#if defined(DNALEN)/* 2-bit packing can be used for sequences only so the define will always be true if pack2bit is used*/
	if (len > DNALEN) return ERR_RUNTIME;
	#endif
	for(i = 0; i < len; ++i) {
		c = seq_nt4_table[(unsigned char)seq[i]];
		if (c < 4) {
			out[i/4] |= c << (2*(3-(i%4)));
		} else {
			return ERR_VALUE;
		}
	}
	return NO_ERROR;
}

unsigned char pDNA8(const char *nib, int len)/*FIXME use AVX instruction set to parallelize this operation*/
{
	int i;
	unsigned char c, v;
	assert(nib != NULL);
	v = 0x00;
	for(i = 0; i < len; ++i) {
		c = seq_nt4_table[(unsigned char)nib[i]];
		if(c < 4) v = (v >> 2) | (c << 6);
	}
	return v;
}
