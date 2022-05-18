#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

#include "sample_main.h"
#include "ketopt.h"
#include "murmur3.h"

void print_sample_help();

enum Error sample_main(int argc, char *argv[]) {
    enum Error err;
    ketopt_t opt;
    int c;
    long int parsed;
    unsigned int i;
    gzFile fp;
    FILE* oh;
    unsigned short r;
    char sep;
    unsigned int seed;
    char buffer[BUFSIZ];
    uint64_t hval[2];

    fp = NULL;
    oh = NULL;
    r = 0;
    sep = '\n';
    seed = 0;

    opt = KETOPT_INIT;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:o:r:p:s:h", longopts)) >= 0) {
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
        } else if (c == 'r') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned short)-1) {
                fprintf(stderr, "Unable to parse sampling rate\n");
                return ERR_OUTOFBOUNDS;
            }
            r = (unsigned short)parsed;
        } else if (c == 'p') {
            if (strlen(opt.arg) != 1) {
                fprintf(stderr, "Separator must be a single character\n");
                return ERR_VALUE;
            }
            sep = opt.arg[0];
        } else if (c == 's') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned int)-1) {
                fprintf(stderr, "Unable to parse option %c\n", c);
                return ERR_OUTOFBOUNDS;
            }
            seed = (unsigned int)parsed;
        } else if (c == 'h') {
            print_sample_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if (r == 0) {
        fprintf(stderr, "Unspecified sampling rate\n");
        return ERR_OPTION;
    }
    if(fp == NULL) {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            fprintf(stderr, "Unable to use stdin as input\n");
            return ERR_OPTION;
        }
    }
    if(oh == NULL) oh = stdout;
    err = NO_ERROR;
    i = 0;
    while((c = gzgetc(fp)) != EOF && err == NO_ERROR) {
        if (c == sep) {
            MurmurHash3_x86_128(buffer, i, seed, &hval);
            if((hval[0] % r) == 0) fprintf(oh, "%.*s\n", i, buffer);
            i = 0;
        } else {
            buffer[i++] = c; /*Ignore k-mers with non-genomic bases (ibf_insert_seq default behaviour)*/
        }        
    }
    if (fp) gzclose(fp);
    return err;
}

void print_sample_help() {
    fprintf(stderr, "[sample] options:\n");
    fprintf(stderr, "\t-i\tinput [stdin]\n");
    fprintf(stderr, "\t-o\toutput [stdout].\n");
    fprintf(stderr, "\t-r\tsampling rate\n");
    fprintf(stderr, "\t-p\tseparator to recognise different input strings [\\n]\n");
    fprintf(stderr, "\t-s\trandom seed [0]\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
