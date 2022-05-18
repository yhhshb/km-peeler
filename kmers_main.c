#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

#include "kmers_main.h"
#include "ketopt.h"
#include "kseq.h"
#include "constants.h"

KSEQ_INIT(gzFile, gzread)

void print_kmers_help();

enum Error kmers_main(int argc, char *argv[]) {
    enum Error err;
    ketopt_t opt;
    int c;
    kseq_t *seq;
    unsigned int i, j;
    gzFile fp;
    FILE* oh;
    long int parsed;
    unsigned char k;

    fp = NULL;
    oh = NULL;
    seq = NULL;
    k = 0;

    opt = KETOPT_INIT;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:o:k:m:S:h", longopts)) >= 0) {
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
            /*silent option that does nothing but is useful to make all comands homegeneous*/
        } else if (c == 'S') {
            /*silent option that does nothing but is useful to make all comands homegeneous*/
        } else if (c == 'h') {
            print_kmers_help();
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
    if(fp == NULL) {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            fprintf(stderr, "Unable to use stdin as input\n");
            return ERR_OPTION;
        }
    }
    if(oh == NULL) {
        oh = stdout;
    }
    err = NO_ERROR;
    seq = kseq_init(fp);
    while(kseq_read(seq) >= 0 && err == NO_ERROR) {
        for(i = j = 0; i < seq->seq.l && err == NO_ERROR; ++i) {
            if (seq_nt4_table[(unsigned char)seq->seq.s[i]] < 4) ++j;
            else j = 0;
            if (j >= k) {
                fprintf(oh, "%.*s\n", k, &seq->seq.s[i-k+1]);
            }
        }
    }
    if (seq) kseq_destroy(seq);
    if (fp) gzclose(fp);
    return err;
}

void print_kmers_help() {
    fprintf(stderr, "[fragment by kmers] options:\n");
    fprintf(stderr, "\t-i\tinput fasta file [stdin]\n");
    fprintf(stderr, "\t-o\toutput file [stdout].\n");
    fprintf(stderr, "\t-k\tk-mer (syncmer) size\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
