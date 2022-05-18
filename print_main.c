#include <stdio.h>
#include "print_main.h"
#include "ketopt.h"
#include "ibflib.h"

void print_print_help();

enum Error print_main(int argc, char** argv) {
    ketopt_t opt;
    ibf_t ibf;
    int c;
    enum Error err;

    opt = KETOPT_INIT;
    err = NO_ERROR;

    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:h", longopts)) >= 0) {
        if (c == 'i') {
            if ((err = ibf_sketch_load(opt.arg, &ibf)) != NO_ERROR) {
                fprintf(stderr, "Unable to read the first invertible bloom filter\n");
                return err;
            }
        } else if (c == 'h') {
            print_print_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    ibf_sketch_print(&ibf, stdout);
    return NO_ERROR;
}

void print_print_help() {
    fprintf(stderr, "[print] options:\n");
    fprintf(stderr, "\t-i\tsketch whose bucket counts are to be printed\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
