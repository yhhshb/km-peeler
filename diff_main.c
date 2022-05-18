#include <stdlib.h>
#include <stdio.h>
#include "diff_main.h"
#include "ketopt.h"
#include "ibflib.h"

void print_diff_help();

enum Error diff_main(int argc, char** argv) {
    ketopt_t opt;
    ibf_t ibf1, ibf2, res;
    int c;
    char *output_path;
    enum Error err;

    opt = KETOPT_INIT;
    output_path = NULL;
    res.seres = NULL;
    res.data = NULL;
    err = NO_ERROR;

    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    while((c = ketopt(&opt, argc, argv, 1, "i:j:o:h", longopts)) >= 0) {
        if (c == 'i') {
            if ((err = ibf_sketch_load(opt.arg, &ibf1)) != NO_ERROR) {
                fprintf(stderr, "Unable to read the first invertible bloom filter\n");
                return err;
            }
        } else if (c == 'j') {
            if ((err = ibf_sketch_load(opt.arg, &ibf2)) != NO_ERROR) {
                fprintf(stderr, "Unable to read the second invertible bloom filter\n");
                return err;
            }
        } else if (c == 'o') {
            output_path = opt.arg;
        } else if (c == 'h') {
            print_diff_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if (output_path == NULL) {
        err = ERR_VALUE;
    }
    if (!err && (err = ibf_sketch_diff(&ibf1, &ibf2, &res)) != NO_ERROR) fprintf(stderr, "Error while computing the ibf difference\n");
    if (!err && (err = ibf_sketch_store(output_path, &res)) != NO_ERROR) fprintf(stderr, "Error saving the difference\n");
    if (!err) err = ibf_sketch_destroy(&ibf1);
    if (!err) err = ibf_sketch_destroy(&ibf2);
    if (!err) err = ibf_sketch_destroy(&res);
    return err;
}

void print_diff_help() {
    fprintf(stderr, "[diff] options:\n");
    fprintf(stderr, "\t-i\tfirst sketch\n");
    fprintf(stderr, "\t-j\tsecond sketch\n");
    fprintf(stderr, "\t-o\tresulting sketch when making $i - $j\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}
