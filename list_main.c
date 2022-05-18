#include <stdio.h>
#include "list_main.h"
#include "ketopt.h"
#include "ibflib.h"

#include <assert.h>

/*If GLEN is defined, then the length of the keys is stored inside the ibf globally*/
#ifdef GLEN /*Nasty trick to bind the global length as a parameter of the print functions declared here*/
keysum_len_t binded_len;
#endif

void print_list_help();

void print_whole_bucket(const bucket_t *bucket, char source, void *unused);
void print_exact_bucket(const bucket_t *bucket, char source, void *unused);

enum Error list_main(int argc, char *argv[]) {
    ketopt_t opt;
    ibf_t ibf;
    unsigned char sketch_read;
    int c;
    enum Error err;

    opt = KETOPT_INIT;
    err = NO_ERROR;
    sketch_read = FALSE;

    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:h", longopts)) >= 0) {
        if (c == 'i') {
            sketch_read = TRUE;
            if ((err = ibf_sketch_load(opt.arg, &ibf)) != NO_ERROR) {
                fprintf(stderr, "Unable to read the first invertible bloom filter\n");
                return err;
            }
        } else if (c == 'h') {
            print_list_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if(!sketch_read) {
        print_list_help();
        return ERR_OPTION;
    }
    #ifdef GLEN
    binded_len = ibf.key_len;
    #endif
    if ((err = ibf_list_seq(&ibf, &print_exact_bucket, NULL), NULL) != NO_ERROR) fprintf(stderr, "Error while peeling the sketch\n");
    ibf_sketch_destroy(&ibf);
    return err;
}

void print_list_help() {
    fprintf(stderr, "[list] options:\n");
    fprintf(stderr, "\t-i\tthe sketch to be listed\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}

void print_whole_bucket(const bucket_t *bucket, char source, void *unused) {
    char sbuf[5];
    unsigned char pack;
    int i, j;
    assert(bucket != NULL);
    sbuf[4] = '\0';
    fprintf(stdout, "%c,", source);
    for (i = 0; i < WSIZE; ++i) {
        pack = bucket->keysum[i];
        for (j = 3; j >= 0; --j) {
            sbuf[j] = seq_nt4_inv_table[pack & INV_MASK];
            pack >>= 2;
        }
        fprintf(stdout, "%s", sbuf);
    }
    #ifndef RPOS
    fprintf(stdout, ",%" format_pos, bucket->position);
    #endif
    fprintf(stdout, "\n");
}

void print_exact_bucket(const bucket_t *bucket, char source, void *unused) {
    keysum_len_t len;
    char sbuf[5];
    unsigned char pack;
    int i, j;
    assert(bucket != NULL);
    #ifndef GLEN
    len = bucket->key_len;
    #elif defined(GLEN)
    len = binded_len;
    #endif
    sbuf[4] = '\0';
    i = 0;
    fprintf(stdout, "%c,", source);
    while (i < len) {
        pack = bucket->keysum[i/4];
        for (j = 3; j >= 0; --j) {
            sbuf[j] = seq_nt4_inv_table[pack & INV_MASK];
            pack >>= 2;
            ++i;
        }
        /*fprintf(stderr, "%d, %u\n", i, bucket->key_len);*/
        if(i > len) {
            sbuf[4 - (i - len)] = '\0';
        }
        fprintf(stdout, "%s", sbuf);
    }
    #ifndef RPOS
    fprintf(stdout, ",%" format_pos, bucket->position);
    #endif
    fprintf(stdout, "\n");
}
