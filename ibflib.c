#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ibflib.h"
#include "murmur3.h"
#include "err.h"
#include "endian_fixer.h"

#include <assert.h>

/*From paper "Invertible Bloom Lookup Tables" (Michael T. Goodrich, Michael Mitzenmacher)*/

#define RMAX 8
enum Access_t {INSERTION, DELETION};
static float ck_table[RMAX] = {0, 0, 0, 1.222, 1.295, 1.425, 1.570, 1.721};

typedef union {
	float f;
	uint32_t u32;
} ufloat32_t;

int ibf_hash_init(unsigned int root_seed, unsigned char r, hash_gen_t **const gen) {
	int i, j;
	assert(gen != NULL);
    if ((*gen = (hash_gen_t*)malloc(r * sizeof(hash_gen_t))) == NULL) return ERR_ALLOC;
    for(i = 0; i < r; ++i) {/*seeds initialization*/
        j = root_seed + i;/*FIXME use a better method?*/
    	(*gen)[i].seed = *( (uint32_t*)((void*)(&j)) );
    }
	return NO_ERROR;
}

int ibf_hash_destroy(hash_gen_t *const gen) {
	assert(gen != NULL);
	if (gen) free(gen);
	return NO_ERROR;
}

int ibf_hash_store(FILE* const out, hash_gen_t const *const gen) {
	/*uint8_t i;*/
	uint32_t buffer32;
	/*uint64_t buffer64;*/
	assert(out != NULL);
	assert(gen != NULL);
	buffer32 = hton32(gen->seed);
	if (fwrite(&buffer32, sizeof buffer32, 1, out) != 1) return ERR_IO;
	return NO_ERROR;/*
	buffer64 = hton64(gen->hash.ls64b);
	fwrite(&buffer64, sizeof buffer64, 1, out);
	buffer64 = hton64(gen->hash.ms64b);
	fwrite(&buffer64, sizeof buffer64, 1, out);*/
}

int ibf_hash_load(FILE *const in, hash_gen_t *const gen) {
	assert(in != NULL);
	assert(gen != NULL);
	if (fread(&gen->seed, sizeof gen->seed, 1, in) != 1) return ERR_IO;
	gen->seed = ntoh32(gen->seed);
	return NO_ERROR;
}

int ibf_hash_seq(char const *const seq, unsigned int start, unsigned int end, unsigned char r, hash_gen_t *const res) {
	uint8_t i;
	assert(seq != NULL);
	assert(res != NULL);
	for(i = 0; i < r; ++i) MurmurHash3_x64_128(&seq[start], end - start, res[i].seed, (void*)( &(res[i].hash) ));
	return NO_ERROR;
}

int ibf_bucket_store(FILE *const out, bucket_t const *const bucket) {
	uint64_t buffer64;
	#ifndef GLEN
	keysum_len_t buffer_len;
	#endif
	#ifndef RPOS
	position_t buffer_pos;
	#endif
	assert(out != NULL);
	assert(bucket != NULL);
	buffer64 = hton64(bucket->counter);
	if (fwrite(&buffer64, sizeof buffer64, 1, out) != 1) return ERR_IO;
	if (fwrite(bucket->keysum, sizeof(uint8_t), WSIZE, out) != WSIZE) return ERR_IO;
	#ifndef GLEN
	buffer_len = hton_len(bucket->key_len);
	if (fwrite(&buffer_len, sizeof(buffer_len), 1, out) != sizeof(buffer_len)) return ERR_IO;
	#endif
	#ifndef RPOS
	buffer_pos = hton_pos(bucket->position);
	if (fwrite(&buffer_pos, sizeof buffer_pos, 1, out) != 1) return ERR_IO;
	#endif
	return NO_ERROR;
}

int ibf_bucket_load(FILE *const in, bucket_t *const bucket) {
	assert(in != NULL);
	assert(bucket != NULL);
	if (fread(&bucket->counter, sizeof bucket->counter, 1, in) != 1) return ERR_IO;
	bucket->counter = ntoh64(bucket->counter);
	if (fread(&bucket->keysum, sizeof(uint8_t), WSIZE, in) != WSIZE) return ERR_IO;
	#ifndef GLEN
	if (fread(&bucket->key_len, sizeof bucket->key_len, 1, in) != 1) return ERR_IO;
	bucket->key_len = ntoh_len(bucket->key_len);
	#endif
	#ifndef RPOS
	if (fread(&bucket->position, sizeof bucket->position, 1, in) != 1) return ERR_IO;
	bucket->position = ntoh_pos(bucket->position);
	#endif
	return NO_ERROR;
}

/*
public access -------------------------------------------------------------------------------------------------
*/

int ibf_sketch_init(unsigned int root_seed, unsigned char r, float epsilon, unsigned int n, ibf_t *const sketch) {
	int err;
	assert(sketch != NULL);
	sketch->repetitions = r;
	sketch->epsilon = epsilon;
	sketch->chunk_size = (unsigned long)ceil((ck_table[sketch->repetitions] + sketch->epsilon) * n / r + 1);
	if ((sketch->data = (bucket_t*)calloc(sketch->chunk_size * sketch->repetitions, sizeof(bucket_t))) == NULL) return ERR_ALLOC;
	if ((err = ibf_hash_init(root_seed, sketch->repetitions, &sketch->seres)) != NO_ERROR) return err;
	return NO_ERROR;
}

int ibf_sketch_copy(ibf_t const *const source, ibf_t *const dest) {
	void* dummy = NULL;
	assert(source != NULL);
	dest->repetitions = source->repetitions;
	dest->epsilon = source->epsilon;
	dest->chunk_size = source->chunk_size;
	if ((dummy = realloc(dest->data, dest->chunk_size * dest->repetitions * sizeof(bucket_t))) != NULL) {
		dest->data = (bucket_t*)dummy;
	} else {
		return ERR_ALLOC;
	}
	if ((dummy = realloc(dest->seres, dest->repetitions * sizeof(hash_gen_t))) != NULL) {
		dest->seres = (hash_gen_t*)dummy;
	} else {
		return ERR_ALLOC;
	}
	if ((dummy = memcpy(dest->data, source->data, dest->chunk_size * dest->repetitions * sizeof(bucket_t))) != dest->data) return ERR_RUNTIME;
	if ((dummy = memcpy(dest->seres, source->seres, dest->repetitions * sizeof(hash_gen_t))) != dest->seres) return ERR_RUNTIME;
	return NO_ERROR;
}

int ibf_sketch_destroy(ibf_t *const sketch) {
	int err;
	assert(sketch != NULL);
	err = NO_ERROR;
	if (sketch->seres != NULL) err = ibf_hash_destroy(sketch->seres);
	if (err == NO_ERROR) sketch->seres = NULL;
	if (sketch->data != NULL) free(sketch->data);
	if (err == NO_ERROR) sketch->data = NULL;
	return err;
}

int ibf_sketch_store(char const *const path, ibf_t const *const sketch) {
	uint64_t i;
	FILE* out;
	ufloat32_t buffer32;
	uint64_t buffer64;
	#ifdef GLEN
	keysum_len_t buffer_len;
	#endif
	assert(path != NULL);
	assert(sketch != NULL);
	if((out = fopen(path, "wb")) == NULL) return ERR_IO;
	if (fwrite(&sketch->repetitions, sizeof sketch->repetitions, 1, out) != 1) return ERR_IO;
	buffer32.f = sketch->epsilon;
	/*buffer32 = *(uint32_t*)(&sketch->epsilon);*/
	buffer32.u32 = hton32(buffer32.u32);
	if (fwrite(&buffer32.u32, sizeof buffer32.u32, 1, out) != 1) return ERR_IO;
	buffer64 = hton64(sketch->chunk_size);
	if (fwrite(&buffer64, sizeof buffer64, 1, out) != 1) return ERR_IO;
	#ifdef GLEN
    buffer_len = hton_len(sketch->key_len);
	if (fwrite(&buffer_len, sizeof buffer_len, 1, out) != 1) return ERR_IO;
    #endif
	
	for(i = 0; i < sketch->chunk_size * sketch->repetitions; ++i) {
		if (ibf_bucket_store(out, &sketch->data[i]) != NO_ERROR) return ERR_IO;
	}
	for(i = 0; i < sketch->repetitions; ++i) {
		if (ibf_hash_store(out, &sketch->seres[i]) != NO_ERROR) return ERR_IO;
	}
	fclose(out);
	return NO_ERROR;
}

int ibf_sketch_load(char const *const path, ibf_t *const sketch) {
	uint64_t i;
	FILE* in;
	ufloat32_t buffer32;
	assert(path != NULL);
	assert(sketch != NULL);
	if ((in = fopen(path, "rb")) == NULL) return ERR_IO;
	if (fread(&sketch->repetitions, sizeof sketch->repetitions, 1, in) != 1) return ERR_IO;
	if (fread(&buffer32.u32, sizeof buffer32.u32, 1, in) != 1) return ERR_IO;
	buffer32.u32 = ntoh32(buffer32.u32);
	/*sketch->epsilon = *(float*)(&buffer32);*/
	sketch->epsilon = buffer32.f;
	if (fread(&sketch->chunk_size, sizeof sketch->chunk_size, 1, in) != 1) return ERR_IO;
	sketch->chunk_size = ntoh64(sketch->chunk_size);
	#ifdef GLEN
	if (fread(&sketch->key_len, sizeof sketch->key_len, 1, in) != 1) return ERR_IO;
	sketch->key_len = ntoh_len(sketch->key_len);
	#endif
	if ((sketch->data = (bucket_t*)malloc(sketch->chunk_size * sketch->repetitions * sizeof(bucket_t))) == NULL) return ERR_ALLOC;
	for(i = 0; i < sketch->chunk_size * sketch->repetitions; ++i) {
		if (ibf_bucket_load(in, &sketch->data[i]) != NO_ERROR) return ERR_IO;
	}
	if ((sketch->seres = (hash_gen_t*)malloc(sketch->repetitions * sizeof(hash_gen_t))) == NULL) return ERR_ALLOC;
	for(i = 0; i < sketch->repetitions; ++i) {
		if (ibf_hash_load(in, &sketch->seres[i]) != NO_ERROR) return ERR_IO;
		sketch->seres[i].hash.ls64b = sketch->seres[i].hash.ms64b = 0;
	}
	fclose(in);
	return NO_ERROR;
}

int ibf_sketch_print(ibf_t const *const sketch, FILE *const strm) {
	uint64_t i, j;
	assert(sketch != NULL);
	for(i = 0; i < sketch->repetitions; ++i) {
		fprintf(strm, "%llu = [", (unsigned long long)i);
		for(j = 0; j < sketch->chunk_size; ++j) {
			if (j != sketch->chunk_size - 1) fprintf(strm, "(%lld"
			#ifndef GLEN
			", %" format_len 
			#endif
			"), ", (unsigned long long)sketch->data[i * sketch->chunk_size + j].counter
			#ifndef GLEN
			, sketch->data[i * sketch->chunk_size + j].key_len
			#endif
			);
			else fprintf(strm, "(%lld"
			#ifndef GLEN
			", %" format_len 
			#endif
			")", (unsigned long long)sketch->data[i * sketch->chunk_size + j].counter
			#ifndef GLEN
			, sketch->data[i * sketch->chunk_size + j].key_len
			#endif
			);
		}
		fprintf(strm, "]\n");
	}
	return NO_ERROR;
}

int ibf_sketch_dump(ibf_t const *const sketch, FILE *const strm) {
	uint64_t i, j, h;
	bucket_t* dumped;
	assert(sketch != NULL);
	fprintf(strm, "r = %u, c = %llu\n", sketch->repetitions, sketch->chunk_size);
	for(i = 0; i < sketch->repetitions; ++i) {
		for(j = 0; j < sketch->chunk_size; ++j) {
			dumped = &sketch->data[i * sketch->chunk_size + j];
			fprintf(strm, "data[%llu, %llu] -> %lld|"
			#ifndef GLEN
			"%" format_len 
			#endif
			"|", i, j, dumped->counter 
			#ifndef GLEN
			,dumped->key_len
			#endif
			);
			for(h = 0; h < WSIZE; ++h) {
				fprintf(strm, "%02X", dumped->keysum[h]);
			}
			fprintf(strm, "\n");
		}
		fprintf(strm, "\n");
	}
	return NO_ERROR;
}

int ibf_sketch_diff(ibf_t const *const a, ibf_t const *const b, ibf_t *const result) {
	int err;
	uint8_t i8;
	unsigned char compatibles;
	uint64_t i, j;
	assert(a != NULL);
	assert(b != NULL);
	assert(result != NULL);
	compatibles = TRUE;
	compatibles &= a->repetitions == b->repetitions;
	compatibles &= a->chunk_size == b->chunk_size;
	for(i8 = 0; i8 < a->repetitions; ++i8) compatibles &= a->seres[i8].seed == b->seres[i8].seed;
	if (!compatibles) {
		return ERR_INCOMPATIBLE;
	}
	if ((err = ibf_sketch_destroy(result)) != NO_ERROR) return err;
	memcpy(result, a, sizeof(ibf_t));
	if ((result->seres = (hash_gen_t*)malloc(a->repetitions * sizeof(hash_gen_t))) == NULL) return ERR_ALLOC;
	memcpy(result->seres, a->seres, a->repetitions * sizeof(hash_gen_t));
	if ((result->data = (bucket_t*)malloc(a->repetitions * a->chunk_size * sizeof(bucket_t))) == NULL) return ERR_ALLOC;
	for(i = 0; i < a->chunk_size * a->repetitions; ++i) {
		result->data[i].counter = a->data[i].counter - b->data[i].counter;
		for(j = 0; j < WSIZE; ++j) result->data[i].keysum[j] = a->data[i].keysum[j] ^ b->data[i].keysum[j];
		#ifndef GLEN
		result->data[i].key_len = a->data[i].key_len ^ b->data[i].key_len;
		#endif
		#ifndef RPOS
		result->data[i].position = a->data[i].position ^ b->data[i].position;
		#endif
	}
	return NO_ERROR;
}

int ibf_access_seq(void const *const seq, unsigned int start, unsigned int end, ibf_t *const sketch, uint8_t *const buffer, uint64_t buffer_len, enum Access_t atype) {
	int i, j, err;
	uint64_t pos;
	assert(seq != NULL);
	assert(start <= end);
	assert(sketch != NULL);
	assert(buffer != NULL);
	if (buffer_len < WSIZE) return ERR_OUTOFBOUNDS;
	memset(buffer, 0, buffer_len);

#if defined(STORE_SEQUENCES) || defined(STORE_VLSEQUENCES) || defined(STORE_FRAGMENTS)/*if buckets contain sequences, init buffer and 2 pack the seq fragment*/
	if (pack2bit(&((char*)seq)[start], end-start, buffer) == NO_ERROR) {/*skip fragments with a base not in {A,C,G,T}*/
#elif defined(STORE_HASHES) /*otherwise, directly copy the hash value (our sequence) into the buffer FIXME, not true anymore*/
	memcpy(buffer, seq, end-start);
	{
#endif
		if((err = ibf_hash_seq((char*)buffer, 0, WSIZE, sketch->repetitions, sketch->seres)) != NO_ERROR) return err;/*hash 2bit sequence*/
		for (j = 0; j < sketch->repetitions; ++j) {
			pos = sketch->seres[j].hash.ls64b % sketch->chunk_size + j * sketch->chunk_size; /*FIXME: lsb is 64 bit, msb is not used to find the position*/
#ifdef DEBUG
			fprintf(stderr, "%.*s,", end - start, &((char*)seq)[start]);
			for(i = 0; i < WSIZE; ++i) fprintf(stderr, "%02X", buffer[i]);
			fprintf(stderr, ",%d,%u,%d,%llu\n", end - start, start, j, sketch->seres[j].hash.ls64b % sketch->chunk_size);
#endif
			switch (atype) {
				case INSERTION:
					++sketch->data[pos].counter;
					break;
				case DELETION:
					--sketch->data[pos].counter;
					break;
				default:
					return ERR_VALUE;
			}
			for(i = 0; i < WSIZE; ++i) {/*add (remove) 2bit-encoded fragment to keysum by XORing*/
				sketch->data[pos].keysum[i] ^= buffer[i];
			}
			#ifndef GLEN
			sketch->data[pos].key_len ^= (keysum_len_t)(end - start); /*FIXME [possible bug] XOR does not work here because many equal lengths by chance*/
			#endif
			#ifndef RPOS
			sketch->data[pos].position ^= (position_t)(start);/*positions are all different*/
			#endif
		}
	}
	return NO_ERROR;
}

int ibf_insert_seq(void const *const seq, int start, int end, ibf_t *const sketch, uint8_t *const buffer, uint64_t buffer_len) {
	return ibf_access_seq(seq, start, end, sketch, buffer, buffer_len, INSERTION);
}

int ibf_delete_seq(void const *const seq, int start, int end, ibf_t *const sketch, uint8_t *const buffer, uint64_t buffer_len) {
	return ibf_access_seq(seq, start, end, sketch, buffer, buffer_len, DELETION);
}

/*
uint64_t find_peelable_bucket(const bucket_t *buckets, uint64_t blen, unsigned char *seen) {
	uint64_t i, c;
	for(i = 0; i < blen ; ++i) {
		if (!seen[i] && ((c = buckets[i].counter) == 1 || c == -1)) return i;
	}
	return blen;
}
*/

unsigned char find_peelable_bucket(bucket_t const *const buckets, uint64_t blen, uint64_t *const last) {
	uint64_t start;
	unsigned char empty;
#ifdef DEBUG
	assert(buckets != NULL);
	assert(last != NULL);
#endif
	start = *last;
	empty = TRUE;
	for(;*last < blen ; ++(*last)) {
		if (buckets[*last].counter == 1 || buckets[*last].counter == -1) return TRUE;
		else if (buckets[*last].counter != 0) empty = FALSE;
	}
	for(*last = 0; *last < start; ++(*last)) {
		if (buckets[*last].counter == 1 || buckets[*last].counter == -1) return TRUE;
		else if (buckets[*last].counter != 0) empty = FALSE;
	}
	*last = blen;
	return empty;
}

#define MAXPASSES 10/*FIXME: can it cause probems?*/

int ibf_list_seq(ibf_t *const sketch, void (*output_bucket)(bucket_t const *const, char, void*), void *iostruct) {
	int err;
	uint64_t i, j, blen, idx, pos, seen;
	unsigned char peelable, too_small;
	assert(sketch != NULL);
	assert(output_bucket != NULL);
	blen = sketch->chunk_size * sketch->repetitions;
	idx = seen = 0;
	too_small = FALSE;
	peelable = find_peelable_bucket(sketch->data, blen, &idx);
	while(idx != blen && seen < MAXPASSES * blen) {/*removed condition !too_small. It now ignores the bucket and continues at idx+1*/
#ifdef DEBUG
		fprintf(stderr, "\n");
		ibf_sketch_print(sketch, stderr);
		fprintf(stderr, "bucket %llu with counter = %lld\n", idx, sketch->data[idx].counter);
#endif
		if ((err = ibf_hash_seq((char*)sketch->data[idx].keysum, 0, WSIZE, sketch->repetitions, sketch->seres)) != NO_ERROR) return err;/*hash 2bit sequence*/
		too_small = TRUE;
		for (j = 0; j < sketch->repetitions; ++j) {/*FIXME: first compute all positions*/
			pos = sketch->seres[j].hash.ls64b % sketch->chunk_size + j * sketch->chunk_size;/*msb not used*/
			if (pos == idx) too_small = FALSE;/*check if idx is in there, if not, skip this index (increment by one)*/
			sketch->seres[j].hash.ls64b = pos;
		}
		if (!too_small) {/*hopefully we will find another bucket which was not an error*/
			for (j = 0; j < sketch->repetitions; ++j) {
				pos = sketch->seres[j].hash.ls64b;
				if (pos != idx) {/*peel all the buckets associated to the found key*/
					sketch->data[pos].counter -= sketch->data[idx].counter;/*update counter*/
					for(i = 0; i < WSIZE; ++i) {/*remove 2bit-encoded fragment to keysum by XORing*/
						sketch->data[pos].keysum[i] ^= sketch->data[idx].keysum[i];
					}
					#ifndef GLEN
					sketch->data[pos].key_len ^= sketch->data[idx].key_len;
					#endif
					#ifndef RPOS
					sketch->data[pos].position ^= sketch->data[idx].position;
					#endif
				}
			}
		} else {
			++idx;
			idx %= blen;
		}
		if (!too_small) {/*now remove the found key itself*/
			output_bucket(&sketch->data[idx], sketch->data[idx].counter == 1 ? 'i' : 'j', iostruct);/*and print it*/
			sketch->data[idx].counter -= sketch->data[idx].counter;/*clear counter of peeled bucket*/
			#ifndef GLEN
			sketch->data[idx].key_len ^= sketch->data[idx].key_len;
			#endif
			#ifndef RPOS
			sketch->data[idx].position ^= sketch->data[idx].position;
			#endif
		}
		++seen;
		peelable = find_peelable_bucket(sketch->data, blen, &idx);
	}
	/*if(too_small) return ERR_OUTOFBOUNDS;*/
	if (!peelable || seen >= MAXPASSES * blen) 
	{
		/*fprintf(stderr, "Warning: unpeelable sketch\n");*/
		return ERR_VALUE;
	}
	return NO_ERROR;
}

int ibf_count_seq(ibf_t const *const sketch, unsigned long *const count) {
	size_t i;
	assert(sketch != NULL);
	assert(count != NULL);
	*count = 0;
	if (sketch->repetitions != 0 || sketch->chunk_size != 0) {
		for(i = 0; i < sketch->chunk_size; ++i) {
			*count += sketch->data[i].counter;
		}
	}
	return NO_ERROR;
}

int ibf_buffer_init(uint8_t **const buffer) {
	assert(buffer != NULL);
	if (*buffer != NULL) {
		free(*buffer);
		*buffer = NULL;
	}
	if ((*buffer = (uint8_t*)malloc(WSIZE)) == NULL) return ERR_ALLOC;
	return NO_ERROR;
}

int ibf_buffer_destroy(uint8_t **const buffer) {
	assert(buffer != NULL);
	if (*buffer) {
		free(*buffer);
		*buffer = NULL;
	}
	return NO_ERROR;
}
