#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <zlib.h>

#include "syncmers_main.h"
#include "ketopt.h"
#include "kseq.h"
#include "kvec2.h"
#include "constants.h"
#include "mmlib.h"

#include <assert.h>

KSEQ_INIT(gzFile, gzread)

void print_syncmers_help();

unsigned int is_fully_genomic(char seq[], int l) {
    unsigned int i;
    for(i = 0; i < l; ++i) {
        if (seq_nt4_table[(uint8_t)seq[i]] > 3) return i;
    }
    return l;
}

enum Error syncmers_main(int argc, char *argv[]) {
    gzFile fp;
    FILE* oh, *dh;
    ketopt_t opt;
    kseq_t *seq;
    enum Error err;
    int c;
    unsigned int i, j, mflen, delta;/*, v;*/
    long int parsed;//, li;
    unsigned char k, m;
    unsigned char segmentation;
    /*char* buffer;*/
    uint32_v_t syncpos;
    uint64_t seed;
    unsigned short grouplen;
    unsigned int extlen;

    fp = NULL;
    dh = NULL;
    oh = NULL;
    seq = NULL;
    k = 0;
    m = 0;
    segmentation = 255;
    /*buffer = NULL;*/
    seed = 42;
    grouplen = 0;
    /*buflen = 0;*/

    opt = KETOPT_INIT;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:o:k:m:g:S:d:h", longopts)) >= 0) {
        if (c == 'i') {
            if ((fp = gzopen(opt.arg, "r")) == NULL) {
                fprintf(stderr, "Unable to open the input file %s\n", opt.arg);
                return ERR_FILE;
            }
        } else if (c == 'o') {
            if ((oh = fopen(opt.arg, "w")) == NULL) {
                fprintf(stderr, "Unable to create output file %s\n", opt.arg);
                return ERR_FILE;
            }
        } else if (c == 'k') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned char)-1) {
                fprintf(stderr, "Unable to parse k-mer length\n");
                return ERR_OUTOFBOUNDS;
            }
            k = (unsigned char)parsed;
        } else if (c == 'm') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > CHAR_MAX || parsed < CHAR_MIN) {
                fprintf(stderr, "Unable to parse minimizer length\n");
                return ERR_OUTOFBOUNDS;
            }
            if (parsed < 0) {
                parsed = -parsed;
                segmentation = TRUE;
            } else segmentation = FALSE;
            m = (unsigned char)parsed;
        } else if (c == 'g') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned short)-1) {
                fprintf(stderr, "Unable to parse group length\n");
                return ERR_OUTOFBOUNDS;
            }
            grouplen = (unsigned short)parsed;
        } else if (c == 'S') {
            seed = strtoull(opt.arg, NULL, 10);
        } else if (c == 'd') {
            if ((dh = fopen(opt.arg, "w")) == NULL) {
                fprintf(stderr, "Unable to create debug file %s\n", opt.arg);
                return ERR_FILE;
            }
        } else if (c == 'h') {
            print_syncmers_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if (k == 0) {
        fprintf(stderr, "Unspecified k\n");
        return ERR_OPTION;
    }
    if (segmentation == 255 || m == 0) m = k;/*m was not specified, default behaviour = output simple k-mers*/
    if (m == k) segmentation = FALSE;/*disable segmentation if m == k (even for m given by user)*/
    mflen = k;
    extlen = k + 1 * grouplen;//If double-extended syncmers set constant to 2
    if (segmentation && grouplen) {
        fprintf(stderr, "m < 0 (segmentation) does not allow grouping of k-mers\n");
        return ERR_OPTION;
    }
    //buflen = 2*k - m;
    if(fp == NULL) {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            fprintf(stderr, "Unable to use stdin as input\n");
            return ERR_OPTION;
        }
    }
    if(oh == NULL) {
        oh = stdout;
    }
    kv_init(syncpos);
    seq = kseq_init(fp);
    err = NO_ERROR;
    /*if(!segmentation) buffer = (char*)malloc(buflen);*/
    while(kseq_read(seq) >= 0 && err == NO_ERROR) {
        syncpos.n = 0;
        if ((err = sync_get_pos(seq->seq.s, seq->seq.l, k, m, seed, &syncpos)) != NO_ERROR) {
            fprintf(stderr, "Error when computing syncmers positions for record: %.*s\n", (int)seq->name.l, seq->name.s);
        }
        if (dh) {
            fprintf(dh, "%.*s\n", seq->seq.l, seq->seq.s);
            for(i = 0; syncpos.n > 0 && i < syncpos.n - 1; ++i) {
                fprintf(dh, "%u,", syncpos.a[i]);
            }
            if (syncpos.n > 0) fprintf(dh, "%u", syncpos.a[syncpos.n-1]);
            fprintf(dh, "\n");
        }
        //for(i = 0; i < syncpos.n; ++i) fprintf(stderr, "%u\n", syncpos.a[i]);/*DEBUG*/
        for(i = j = 0; i < syncpos.n && err == NO_ERROR; ++i) {
            /*fprintf(stderr, "%u\n", mmpos.a[i]);*/
            if (segmentation) {/*segments*/
                delta = syncpos.a[i] - j;
                if (delta <= k) fprintf(oh, "%.*s\n", delta, &seq->seq.s[j]);
            } else if (grouplen) {/*extended syncmers*/
                if ((syncpos.a[i] + k + grouplen) < seq->seq.l) {// && syncpos.a[i] >= grouplen) {
                    //fprintf(stderr, "%u, %u, pos = %u\n", grouplen, extlen, syncpos.a[i] - grouplen);
                    fprintf(oh, "%.*s\n", extlen, &seq->seq.s[syncpos.a[i]]);// - grouplen]);
                }
            } else {/*simple syncmers*/
                fprintf(oh, "%.*s\n", mflen, &seq->seq.s[syncpos.a[i]]);
            }
            j = syncpos.a[i];
        }
        if (segmentation) { /*Add last segment, if needed*/
            delta = seq->seq.l - j;
            if (delta <= k) fprintf(oh, "%.*s\n", (int)(delta), &seq->seq.s[j]);
        } else if (grouplen) {
            /*
            Quick and dirty solution, one should scan the sequence in order to find N's, see commented blocks at the end of this file.
            Unfortunately, I am in a hurry and I don't have time to do it properly.
            */
            fprintf(oh, "%.*s\n", extlen, seq->seq.s);/*print first extended k-mer*/
            fprintf(oh, "%.*s\n", extlen, &seq->seq.s[seq->seq.l - extlen]);/*print the last extended k-mer*/
        }
    }
    if (seq) kseq_destroy(seq);
    kv_destroy(syncpos);
    /*if (buffer) free(buffer);*/
    gzclose(fp);
    return err;
}

void print_syncmers_help() {
    fprintf(stderr, "[fragment by syncmers] options:\n");
    fprintf(stderr, "\t-i\tinput fasta file [stdin]\n");
    fprintf(stderr, "\t-o\toutput file [stdout].\n");
    fprintf(stderr, "\t-k\tk-mer (syncmer) size\n");
    fprintf(stderr, "\t-m\tminimizer size for finding syncmers [k]. \n\t\tm > 0 output each syncmer.\n\t\tm < 0 fragments the input sequence at syncmer positions\n");
    fprintf(stderr, "\t-g\tif a syncmer is found at position p, print seq[p:p+k+g]. Incompatible with m < 0 [0]\n");
    fprintf(stderr, "\t-S\tseed [42]\n");
    fprintf(stderr, "\t-d\tdebug file\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}

/*
else if (grouplen) {//print blocks at the beginning and end
            j = 0;
            li = 0;
            //print fragment starting from 0 to the first syncmer or interruption
            if (syncpos.n > 0) {
                if (syncpos.a[j] != 0 && is_fully_genomic(&seq->seq.s[0], syncpos.a[j]) == syncpos.a[j]) 
                    fprintf(oh, "%.*s\n", syncpos.a[j], &seq->seq.s[0]);//if there is an inside the first block it is handled by the following loop
                for(i = 0; i < seq->seq.l; ++i) {//re-scan the sequence looking for N's
                    if (seq_nt4_table[(uint8_t)seq->seq.s[i]] > 3)  { //found an interruption
                        //find the positions of the previous and next syncmers
                        for(; i > syncpos.a[j]; ++j) {}
                        //print fragment starting from previous (extended) syncmer to the position of the interruption
                        if (j > 0) fprintf(oh, "%.*s\n", i - (syncpos.a[j - 1] + extlen), &seq->seq.s[syncpos.a[j - 1] + extlen]);
                        //now print the fragment starting from the interruption to the next syncmer
                        if (syncpos.a[j] != i + 1) fprintf(oh, "%.*s\n", syncpos.a[j] - i - 1, &seq->seq.s[i + 1]);
                    }
                }
            }
            for(li = syncpos.n - 1, parsed = 0; li > 0 && (parsed = seq->seq.l - (syncpos.a[li] + extlen)) < 0; --li) {}
            fprintf(stderr, "li = %lu --> syncmer at %u | parsed = %ld\n", li, syncpos.a[li], parsed);
            if (li >= 0) {
                if (parsed) fprintf(oh, "%.*s\n", parsed, &seq->seq.s[syncpos.a[li] + extlen]);
            } else {
                fprintf(oh, "%.*s\n", seq->seq.l, seq->seq.s);
            }
        }
*/

/*
for(li = syncpos.n - 1; li > 0 && (parsed = seq->seq.l - (syncpos.a[li] + extlen)) < 0; --li) {}
if (li >= 0) {
    fprintf(oh, "%.*s", parsed, &seq->seq.s[syncpos.a[li] + extlen]);
    for(i = 0; i < extlen - parsed; ++i) { //pad last fragment with A's
        fprintf(oh, "A");
    }
    fprintf(oh, "\n");
}
*/

/*
if ((v = is_fully_genomic(&seq->seq.s[i], extlen)) == extlen) { print segment just after it
    fprintf(oh, "%.*s\n", extlen, &seq->seq.s[i]);
} else {skip N and print block just before N
    i += v;
    if (is_fully_genomic(&seq->seq.s[i - extlen], extlen) == extlen) fprintf(oh, "%.*s\n", extlen, &seq->seq.s[i - extlen]);
    j = 0;
}
*/

/*fprintf(oh, "%.*s\n", extlen, &seq->seq.s[seq->seq.l - (seq->seq.l > extlen ? extlen : seq->seq.l)]);*/ /*last extended syncmer at the end*/
