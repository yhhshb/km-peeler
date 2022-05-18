#include <stdio.h>
#include "dump_main.h"
#include "ketopt.h"
#include "ibflib.h"

void print_dump_help();

enum Error dump_main(int argc, char** argv) {
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
                fprintf(stderr, "Unable to read the invertible bloom filter\n");
                return err;
            }
        } else if (c == 'h') {
            print_dump_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if(!sketch_read) {
        print_dump_help();
        return ERR_OPTION;
    }
    ibf_sketch_dump(&ibf, stdout);
    return NO_ERROR;
}

void print_dump_help() {
    fprintf(stderr, "[dump] options:\n");
    fprintf(stderr, "\t-i\tsketch whose bucket counts are to be printed\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
