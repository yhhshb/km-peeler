#include "build_main.h"
#include <stdlib.h>
#include <stdio.h>

#include "ketopt.h"
#include "constants.h"
#include "ibflib.h"

#include <assert.h>

int check_build_args(unsigned int n, unsigned char r, float e, char *opath);
void print_build_help();

/*
 * Construction algorithm for an IBF built on a set of k-mers.
 * The input must be in human-readable form.
 * Each k-mer must be at the beginning of a line.
 */
enum Error build_main(int argc, char *argv[]) {
    FILE *fp;
    char *output_path;
    int c, i, err;
    unsigned char r, l;
    unsigned int n, s;
    long parsed;
    float e;
    ketopt_t opt;
    char* kmer;
    uint8_t *ibfbuf;
    uint64_t blen;
    ibf_t ibf;

    assert(argv != NULL);

    opt = KETOPT_INIT;
    fp = NULL;
    output_path = NULL;
    r = 3;
    n = 0;
    s = 42;
    e = 0;
    err = NO_ERROR;
    kmer = NULL;
    ibfbuf = NULL;
    blen = 0;
    l = 0;

    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:o:n:r:e:s:l:h", longopts)) >= 0) {
        if (c == 'i') {
            if ((fp = fopen(opt.arg, "r")) == NULL) {
                fprintf(stderr, "Unable to open the input file\n");
                return ERR_FILE;
            }
        } else if (c == 'o') {
            /*if ((output_path = (char*)malloc(strlen(opt.arg) + 1)) == NULL) {
                fprintf(stderr, "Unable to allocate space for output filename\n");
                return ERR_ALLOC;
            }
            strcpy(output_path, opt.arg);*/
            output_path = opt.arg;
        } else if (c == 'l') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned char)-1) {
                fprintf(stderr, "Unable to parse option %c\n", c);
                return ERR_OUTOFBOUNDS;
            }
            l = (unsigned char)parsed;
        } else if (c == 'n') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned int)-1) {
                fprintf(stderr, "Unable to parse option %c\n", c);
                return ERR_OUTOFBOUNDS;
            }
            n = (unsigned int)parsed;
        } else if (c == 'r') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned char)-1) {
                fprintf(stderr, "Unable to parse option %c\n", c);
                return ERR_OUTOFBOUNDS;
            }
            r = (unsigned char)parsed;
        } else if (c == 'e') {
            e = atof(opt.arg);
        } else if (c == 's') {
            parsed = strtol(opt.arg, NULL, 10);
            if (parsed > (unsigned int)-1) {
                fprintf(stderr, "Unable to parse option %c\n", c);
                return ERR_OUTOFBOUNDS;
            }
            s = (unsigned int)parsed;
        } else if (c == 'h') {
            print_build_help();
            /*if (output_path != NULL) free(output_path);*/
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            /*if (output_path != NULL) free(output_path);*/
            return ERR_OPTION;
        }
    }
    if (!check_build_args(n, r, e, output_path)) return ERR_OPTION;
    if(fp == NULL) fp = stdin;

    if (err == NO_ERROR) {
        err = ibf_buffer_init(&ibfbuf);
        blen = WSIZE;
        print_error(err, "buffer init");
    }
    if (err == NO_ERROR) {
        err = ibf_sketch_init(s, r, e, n, &ibf);
        print_error(err, "sketch init");
    }

    if (err == NO_ERROR) {
        kmer = (char*)malloc(BUFSIZ);
        if (kmer == NULL) {
            err = ERR_ALLOC;
            print_error(err, "kmer buffer allocation");
        }
    }
    if (err == NO_ERROR) {
#ifdef DEBUG
        fprintf(stderr, "seq,2bit,length,position,row,col\n");
#endif
        i = 0;
        #ifdef GLEN
        ibf.key_len = 0;
        #endif
        while((c = fgetc(fp)) != EOF && err == NO_ERROR) {
            switch (c) {
                case '\n':
                    /*if (i != k) err = ERR_RUNTIME;*/
                    if (err == NO_ERROR && i > 0) err = ibf_insert_seq(kmer, 0, i, &ibf, ibfbuf, blen);
                    if (err == NO_ERROR && l != 0) if (i > l) err = ERR_VALUE;
                    #ifdef GLEN
                    if (ibf.key_len < i) ibf.key_len = i;
                    #endif
                    i = 0;
                    break;
                default:
                    kmer[i++] = c; /*Ignore k-mers with non-genomic bases (ibf_insert_seq default behaviour)*/
            }
            
        }
    }

    if (kmer) free(kmer);
    if (ibfbuf) ibf_buffer_destroy(&ibfbuf);
    if (fp) fclose(fp);
    /*ibf_sketch_dump(&ibf, stderr);*/
    if (err == NO_ERROR) {
        err = ibf_sketch_store(output_path, &ibf);
        if (err != NO_ERROR) print_error(err, "IBF save");
    }
    if (err == NO_ERROR) {
        err = ibf_sketch_destroy(&ibf);
        if (err != NO_ERROR) print_error(err, "sketch destroy");
    }
    /*if (output_path) free(output_path);*/
    return err;
}

int check_build_args(unsigned int n, unsigned char r, float e, char *opath) {
    if(n == 0) {
        fprintf(stderr, "Unspecified n\n");
        return FALSE;
    }
    if(r < 3 || r > 7) {
        fprintf(stderr, "2 < r < 8\n");
        return FALSE;
    }
    if(e < 0.0) {
        fprintf(stderr, "0 <= e <= 1\n");
        return FALSE;
    }
    if(opath == NULL) {
        fprintf(stderr, "Output file not specified\n");
        return FALSE;
    }
    return TRUE;
}

void print_build_help() {
    fprintf(stderr, "[build] options:\n");
    fprintf(stderr, "\t-i\tinput set of k-mers [stdin]\n");
    fprintf(stderr, "\t-o\tInvertible Bloom Filter file (binary output)\n");
    fprintf(stderr, "\t-n\tnumber of differences to track (0 < n)\n");
    fprintf(stderr, "\t-r\tnumber of hash functions [3] (3 <= r <= 7)\n");
    fprintf(stderr, "\t-e\tepsilon [0] (0 <= epsilon)\n");
    fprintf(stderr, "\t-s\trandom seed [42]\n");
    fprintf(stderr, "\t-l\tmaximum length of input sequences, used for checking correctness\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
