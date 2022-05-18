#include <stdlib.h>
#include <stdio.h>
#include "jaccard_main.h"
#include "ketopt.h"
#include "ibflib.h"
/*#include "minHash.h"*/

typedef struct {
    unsigned long long unique_i_size;
    unsigned long long unique_j_size;
} callback_t;

void print_jaccard_help();

void increment_difference_size(const bucket_t *bucket, char source, void *sizes);

enum SketchType {IBF, MINHASH};

int compute_ibf_jaccard(char const *const ibf1_path, char const *const ibf2_path);

enum Error jaccard_main(int argc, char** argv) {
    ketopt_t opt;
    int c;
    char *path1, *path2;
    enum Error err;
    enum SketchType stype;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};

    stype = IBF;/* IBF is always used because y option is not available anymore */
    err = NO_ERROR;
    opt = KETOPT_INIT;
    path1 = path2 = NULL;
    while((c = ketopt(&opt, argc, argv, 1, "i:j:h", longopts)) >= 0) {
        if (c == 'i') {
            path1 = opt.arg;
        } else if (c == 'j') {
            path2 = opt.arg;
        } else if (c == 'y') {
            if (strcmp(opt.arg, "ibf") == 0) stype = IBF;
            else if (strcmp(opt.arg, "mh") == 0) stype = MINHASH;
            else return ERR_OPTION;
        } else if (c == 'h') {
            print_jaccard_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if (!path1 || !path2) return ERR_OPTION;
    switch (stype) {
        case IBF:
            if (!err) err = compute_ibf_jaccard(path1, path2);
            break;
        /*
        case MINHASH:
            if (!err) err = compute_minhash_jaccard(path1, path2);
            break;
        */
        default:
            fprintf(stderr, "[Jaccard] This should never happen\n");
    }
    return err;
}

void print_jaccard_help() {
    fprintf(stderr, "[jaccard] options:\n");
    fprintf(stderr, "\t-i\tfirst sketch\n");
    fprintf(stderr, "\t-j\tsecond sketch\n");
    /*fprintf(stderr, "\t-y\tsketch type (ibf, mh) [ibf]\n");*/
    fprintf(stderr, "\t-h\tshow this help\n");
}

void increment_difference_size(const bucket_t *bucket, char source, void *sizes) {
    callback_t *dummy = (callback_t*) sizes;
    if (source == 'i') ++dummy->unique_i_size;
    else ++dummy->unique_j_size;
}

int compute_ibf_jaccard(char const *const ibf1_path, char const *const ibf2_path) {
    int err;
    unsigned long L0i, L0j;
    double jaccard, containment_i_j, containment_j_i;
    ibf_t ibf1, ibf2, res;
    callback_t increments;
    res.seres = NULL;
    res.data = NULL;
    increments.unique_i_size = 0;
    increments.unique_j_size = 0;
    err = NO_ERROR;
    if (!err && (err = ibf_sketch_load(ibf1_path, &ibf1)) != NO_ERROR) fprintf(stderr, "Unable to read the first invertible bloom filter\n");
    if (!err && (err = ibf_sketch_load(ibf2_path, &ibf2)) != NO_ERROR) fprintf(stderr, "Unable to read the second invertible bloom filter\n");
    if (!err) if ((err = ibf_sketch_diff(&ibf1, &ibf2, &res)) != NO_ERROR) fprintf(stderr, "Error while computing the ibf difference\n");
    if (!err) err = ibf_count_seq(&ibf1, &L0i);
    if (!err) err = ibf_count_seq(&ibf2, &L0j);
    if (!err) err = ibf_sketch_destroy(&ibf1);
    if (!err) err = ibf_sketch_destroy(&ibf2);
    
    if (!err) if ((err = ibf_list_seq(&res, &increment_difference_size, &increments)) != NO_ERROR) fprintf(stderr, "Error while peeling the sketch for jaccard computation\n");
    if (!err) err = ibf_sketch_destroy(&res);
    if (!err) jaccard = ((double)(L0i - increments.unique_i_size)) / (L0i + increments.unique_j_size);
    if (!err) if (jaccard != ((double)(L0j - increments.unique_j_size)) / (L0j + increments.unique_i_size)) {
        fprintf(stderr, "Inconsistent Jaccard\n");
        err = ERR_RUNTIME;
    }
    if (!err) {
        containment_i_j = (double)(L0i - increments.unique_i_size) / L0i;
        containment_j_i = (double)(L0j - increments.unique_j_size) / L0j;
        fprintf(stdout, "%f,%f,%f\n", jaccard, containment_i_j, containment_j_i);
    } else fprintf(stdout, "NaN,NaN,NaN\n");
    return err;
}

void print_count_help() {
    fprintf(stderr, "[count] options:\n");
    fprintf(stderr, "\t-i\tibf sketch to count\n");
    fprintf(stderr, "\t-h\tshow this help\n");
}

enum Error count_main(int argc, char** argv) {
    ketopt_t opt;
    int c;
    char *path;
    enum Error err;
    unsigned long L0;
    ibf_t ibf;
    static ko_longopt_t longopts[] = {{NULL, 0, 0}};
    err = NO_ERROR;
    opt = KETOPT_INIT;
    path = NULL;
    while((c = ketopt(&opt, argc, argv, 1, "i:j:h", longopts)) >= 0) {
        if (c == 'i') {
            path = opt.arg;
        } else if (c == 'h') {
            print_count_help();
            return NO_ERROR;
        } else {
            fprintf(stderr, "Option -%c not available\n", c);
            return ERR_OPTION;
        }
    }
    if (!path) return ERR_OPTION;
    L0 = 0;
    if (!err && (err = ibf_sketch_load(path, &ibf)) != NO_ERROR) fprintf(stderr, "Unable to load the invertible bloom filter\n");
    if (!err) err = ibf_count_seq(&ibf, &L0);
    if (!err) err = ibf_sketch_destroy(&ibf);
    if (!err) fprintf(stdout, "%lu", L0);
    return err;
}

/* in-house minhash sketches are not really needed, use Mash instead */
/*
int compute_minhash_jaccard(char const *const ibf1_path, char const *const ibf2_path) {
    int err;
    minhash_t mh1, mh2;
    double jaccard;
    err = NO_ERROR;
    if (!err && (err = minhash_sketch_load(ibf1_path, &mh1)) != NO_ERROR) fprintf(stderr, "Unable to read the first minHash sketch\n");
    if (!err && (err = minhash_sketch_load(ibf2_path, &mh2)) != NO_ERROR) fprintf(stderr, "Unable to read the second minHash sketch\n");
    if (!err && (err = minhash_compare(&mh1, &mh2, &jaccard)) != NO_ERROR) fprintf(stderr, "Error while computing the minHash jaccard\n");
    if (!err) err = minhash_sketch_destroy(&mh1);
    if (!err) err = minhash_sketch_destroy(&mh2);
    if (!err) fprintf(stdout, "%f,NaN,NaN\n", jaccard);
    else fprintf(stdout, "NaN,NaN,NaN\n");
    return err;
}
*/
