#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdint.h>

/********************************************************************************/

/*Target application selection*/

#include "compile_options.h"

#define CEILING(x,y) (((x) + (y) - 1) / (y))

#if defined(STORE_HASHES) && !defined(STORE_SEQUENCES) /*if we do not want to store sequences, store their 64bit hashes instead*/
    #define HASHLEN STORE_HASHES
    #define WSIZE ((CEILING(HASHLEN,8)))
    #define GLEN                 /*Global sequence length information*/
    #define RPOS                 /*remove positions from ibf struct*/
    #define MMML 0
    #define MAXSEQLEN 0          /*No positions in the buckets*/
#elif defined(STORE_SEQUENCES) && !defined(STORE_HASHES)
    #define DNALEN STORE_SEQUENCES
    #define WSIZE ((CEILING(2*DNALEN,8)))
    #define GLEN                 /*Global sequence length information*/
    #define RPOS                 /*remove positions from ibf struct*/
    #define MMML 0
    #define MAXSEQLEN 0          /*No positions in the buckets*/
#elif defined(STORE_FRAGMENTS)
    #define DNALEN STORE_FRAGMENTS
    #define WSIZE ((CEILING(2*DNALEN,8)))
    #define MMML 32              /*Max mm length*/
    #define MAXSEQLEN 0xFFFFFFFF /*Maximum allowed sequence length (store fragment position in buckets)*/
#elif defined(STORE_VLSEQUENCES)
    #define DNALEN STORE_VLSEQUENCES
    #define WSIZE ((CEILING(2*DNALEN,8)))
    #define RPOS                 /*remove positions from ibf struct*/
    #define MMML 0
    #define MAXSEQLEN 0          /*No positions in the buckets*/
#endif

#define MWW (DNALEN - MMML)      /*Maximum Window Width*/

/********************************************************************************/

/*Utilities and constants*/

#define TRUE 1
#define FALSE 0
#define INV_MASK 0x03

#define parse_check(unsigned_type, opt, buffer_var, dealloc, retval) do {\
        if (buffer_var > (unsigned_type)-1) {\
            fprintf(stderr, "Unable to parse option %c\n", opt);\
            if (dealloc) {\
                free(dealloc);\
                dealloc = NULL;\
            }\
            retval = FALSE;\
        } else {\
            retval = TRUE;\
        }\
    } while(0)

unsigned char seq_nt4_table[256];

char seq_nt4_inv_table[5];

int pack2bit(const char *seq, unsigned char len, unsigned char *out);

#endif/*CONSTANTS_H*/
