#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <zlib.h>

#include "minimizers_main.h"
#include "ketopt.h"
#include "kseq.h"
#include "kvec2.h"
#include "constants.h"
#include "mmlib.h"

#include <assert.h>

KSEQ_INIT(gzFile, gzread)

void print_minimizers_help();

enum Error minimizers_main(int argc, char *argv[]) {
    gzFile fp;
    FILE* oh;
    ketopt_t opt;
    kseq_t *seq;
    enum Error err;
    int c;
    unsigned int i, j, mflen, buflen, delta;
    long int parsed;
    unsigned char k, m;
    unsigned short w;
    unsigned char segmentation, split;/*, canonical;*/
    char* buffer;
    uint32_v_t mmpos;
    uint64_t seed;

    fp = NULL;
    oh = NULL;
    seq = NULL;
    k = 0;
    m = 0;
    w = 0;
    segmentation = 255;
    /*canonical = FALSE;*/
    split = FALSE;
    buffer = NULL;
    buflen = 0;
    seed = 42;

    opt = KETOPT_INIT;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:o:k:m:w:S:sh", longopts)) >= 0) {
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
        } else if (c == 's') {
            split = TRUE;
        } else if (c == 'w') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned short)-1) {
                fprintf(stderr, "Unable to parse window length\n");
                return ERR_OUTOFBOUNDS;
            }
            w = (unsigned short)parsed;
        } else if (c == 'S') {
            seed = strtoull(opt.arg, NULL, 10);
        } else if (c == 'c') {
            /*canonical = TRUE;*//*not used*/
        } else if (c == 'h') {
            print_minimizers_help();
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
    if (w) {
        if (m) {/*other error conditions when w is active*/
            fprintf(stderr, "options <m> and <w> are mutually exclusive\n");
            return ERR_OPTION;
        }
        /*special settings when w is active should go here (it is empty as of now)*/
        segmentation = TRUE;/*disables buffer allocation (see later) and split option*/
    }
    if (segmentation == 255 || m == 0) m = k;/*m was not specified, default behaviour = output simple k-mers*/
    if (m == k) segmentation = FALSE;/*disable segmentation if m == k (even for m given by user)*/
    if (split) mflen = k;
    else mflen = 2*k - m;
    buflen = 2*k - m;
    if(fp == NULL) {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            fprintf(stderr, "Unable to use stdin as input\n");
            return ERR_OPTION;
        }
    }
    if(oh == NULL) {
        oh = stdout;
    }

    kv_init(mmpos);
    seq = kseq_init(fp);
    err = NO_ERROR;
    if(!segmentation) buffer = (char*)malloc(buflen);
    while(kseq_read(seq) >= 0 && err == NO_ERROR) {
        mmpos.n = 0;
        if (w) {
            if ((err = mm_get_pos(seq->seq.s, seq->seq.l, k, w, seed, &mmpos)) != NO_ERROR) 
                fprintf(stderr, "Error when computing minimizer positions for record: %.*s\n", (int)seq->name.l, seq->name.s);
            for(i = 0; i < mmpos.n && err == NO_ERROR; ++i) {
                fprintf(oh, "%.*s\n", k, &seq->seq.s[mmpos.a[i]]);
            }
        } else {
            if ((err = mm_get_pos(seq->seq.s, seq->seq.l, m, k-m+1, seed, &mmpos)) != NO_ERROR) 
                fprintf(stderr, "Error when computing minimizer positions for record: %.*s\n", (int)seq->name.l, seq->name.s);
            /*for(i = 0; i < mmpos.n; ++i) fprintf(stderr, "%u\n", mmpos.a[i]);*//*DEBUG*/
            for(i = j = 0; i < mmpos.n && err == NO_ERROR; ++i) {
                /*fprintf(stderr, "%u\n", mmpos.a[i]);*/
                if (segmentation) {
                    delta = mmpos.a[i] - j;
                    if (delta <= k) fprintf(oh, "%.*s\n", delta, &seq->seq.s[j]);
                } else {/*splitted or grouped syncmers*/
                    if (mmpos.a[i] < (k - m)) {/*first mm positions, pad with A's*/
                        memset(buffer, 'A', buflen);/*buffer is always 2*k-m long*/
                        memcpy(&buffer[k - m - mmpos.a[i]], seq->seq.s, k + mmpos.a[i]);
                        if (split) fprintf(oh, "%.*s\n", mflen, buffer + k-m);
                        fprintf(oh, "%.*s\n", mflen, buffer);
                    } else if (seq->seq.l - mmpos.a[i] < k) {/*last positions closer than k to the end*/
                        memset(buffer, 'A', buflen);
                        memcpy(buffer, &seq->seq.s[mmpos.a[i]-k+m], seq->seq.l - mmpos.a[i]+k-m);// - mmpos.a[i]);
                        if (split) fprintf(oh, "%.*s\n", mflen, buffer +k-m);
                        fprintf(oh, "%.*s\n", mflen, buffer);
                    } else {/*padding not needed here*/
                        if (split) fprintf(oh, "%.*s\n", mflen, &seq->seq.s[mmpos.a[i]]);
                        fprintf(oh, "%.*s\n", mflen, &seq->seq.s[mmpos.a[i]-k+m]);
                    }
                }
                j = mmpos.a[i];
            }
            if (segmentation) { /*Add last segment, if needed*/
                delta = seq->seq.l - j;
                if (delta <= k) fprintf(oh, "%.*s\n", (int)(delta), &seq->seq.s[j]);
            }/*
            } else {
                memset(buffer, 'A', 2*k - m);
                fprintf(stderr, "k = %u, m = %u, j = %u, len = %lu\n", (uint32_t)k, (uint32_t)m, j, seq->seq.l);
                memcpy(buffer, &seq->seq.s[j], seq->seq.l - j);
                if (split) fprintf(oh, "%.*s\n", mflen, buffer + m);
                fprintf(oh, "%.*s\n", mflen, buffer);
            }*/
        }
    }
    if (seq) kseq_destroy(seq);
    kv_destroy(mmpos);
    if (buffer) free(buffer);
    return err;
}

void print_minimizers_help() {
    fprintf(stderr, "[fragment by minimizers] options:\n");
    fprintf(stderr, "\t-i\tinput fasta file [stdin]\n");
    fprintf(stderr, "\t-o\toutput file [stdout].\n");
    fprintf(stderr, "\t-k\tk-mer size\n");
    fprintf(stderr, "\t-m\tminimizer size for grouping [k]. \n\t\tm > 0 output [p-(k-m):p+k] bases for each minimizer at position p.\n\t\tm < 0 fragments the input sequence at minimizer positions, disables option <s>\n");
    fprintf(stderr, "\t-s\tsplit each group of k-mers into its constituent syncmers <inactive>\n");
    fprintf(stderr, "\t-w\twindow length for indexing (number of k-mers), disables option <s>, incompatible with option <m>\n");
    fprintf(stderr, "\t-S\tseed [42]\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
