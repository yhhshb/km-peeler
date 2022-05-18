#include <sys/stat.h>
#include <zlib.h>
#include <stdio.h>

#include "kseq.h"
#include "kvec2.h"
#include "constants.h"
#include "minHash.h"
#include "ibflib.h"
#include "collection_main.h"

#include <assert.h>

KSEQ_INIT(gzFile, gzread)

typedef struct {
	size_t n;
	size_t m;
	minhash_t *a;
} minhash_v_t; /* collection minhash of sketches */

typedef char** path_itr;

typedef struct {
    uint32_t opt_end; /* index pointing after the last found option */
    uint8_t k;  /* k-mer length */
    uint8_t w;  /* hash width in number of bits */
    uint32_t s; /* global seed */
    uint64_t z; /* number of hashes to be stored inside each sketch */
    uint64_t n;
    uint8_t r;
    float e;
    char *opath;
    path_itr alice_start, alice_stop, bob_start, bob_stop;
} param_t; /* Here ketopt.h is not used because we want to replicate the behaviour of Mash */

typedef struct {
    minhash_v_t *alice, *bob;
    uint8_t *dup_alice, *dup_bob;
} callback_t;

void print_collection_help();
void init_options(param_t *parameters);
static enum Error option_parser(int *opt_idx, char opt, char **arg, param_t *parameters);
enum Error parse_options(int argc, char *argv[], char const * const opt_string, enum Error (*parser)(int*, char, char**, param_t*), param_t *params);
enum Error check_options(param_t const * const parameters);
enum Error fill_minhash_sketch(param_t const * const opts, path_itr start, path_itr end, minhash_v_t * const mhdb);
enum Error mark_duplicates(minhash_v_t const * const mhdb, uint8_t ** const duplicate);
void get_sketches(bucket_t const * const bucket, char origin, void *callback_io);

enum Error collection_main(int argc, char *argv[]) {
    int i, j, err;
    param_t opts;
    minhash_v_t alice, bob;
    minhash_t *a, *b;
    //minhash_t dummy;
    uint8_t *duplicate_alice, *duplicate_bob;
    uint8_t *ibfbuf;
    uint64_t buflen;
    ibf_t ibf_alice, ibf_bob, diff;
    callback_t callback_io;
    assert(argv);
    err = 0;
    ibfbuf = NULL;
    buflen = 0;
    #ifdef GLEN
    ibf_alice.key_len = ibf_bob.key_len = 0;
    #endif
    //dummy.hash_width = 0;/* suppress warning */
    init_options(&opts);
    // fprintf(stderr, "[Log] Begin\n");
    if (!err) err = parse_options(argc, argv, "k:z:w:n:r:e:s:o:a:b:h", option_parser, &opts);
    if (!err) err = check_options(&opts);
    if (err) return err;

    // for(i = 0; i < opts.alice_stop - opts.alice_start; ++i) fprintf(stderr, "%s\n", *(opts.alice_start + i));
    // for(i = 0; i < opts.bob_stop - opts.bob_start; ++i) fprintf(stderr, "%s\n", *(opts.bob_start + i));

    kv_init(alice);
    kv_init(bob);
    // fprintf(stderr, "index at position %llu after options\n", opts.opt_end);
    // fprintf(stderr, "k = %d, w = %d, s = %u, z = %llu\n", opts.k, opts.w, opts.s, opts.z);
    if (!err) err = fill_minhash_sketch(&opts, opts.alice_start, opts.alice_stop, &alice);
    if (!err) err = fill_minhash_sketch(&opts, opts.bob_start, opts.bob_stop, &bob);
    /* Mark duplicate sketches */
    duplicate_alice = NULL;
    duplicate_bob = NULL;
    if (!err) err = mark_duplicates(&alice, &duplicate_alice);
    if (!err) err = mark_duplicates(&bob, &duplicate_bob);
    // fprintf(stderr, "[Log] duplicates marked\n");
    /* Insert sketches into IBF */
    if (!err) err = ibf_buffer_init(&ibfbuf);
    buflen = WSIZE;
    // fprintf(stderr, "----- WSIZE = %llu\n", buflen);
    if (!err) err = ibf_sketch_init(opts.s, opts.r, opts.e, opts.n, &ibf_alice);
    for(i = 0; !err && i < alice.n; ++i) {
        if (!duplicate_alice[i]) err = ibf_insert_seq(alice.a[i].hashes, 0, alice.a[i].size * alice.a[i].hash_width, &ibf_alice, ibfbuf, buflen);
    }
    if (!err) err = ibf_sketch_init(opts.s, opts.r, opts.e, opts.n, &ibf_bob);
    for(i = 0; !err && i < bob.n; ++i) {
        if (!duplicate_bob[i]) err = ibf_insert_seq(bob.a[i].hashes, 0, bob.a[i].size * bob.a[i].hash_width, &ibf_bob, ibfbuf, buflen);
    }
    if (ibfbuf) ibf_buffer_destroy(&ibfbuf);
    // fprintf(stderr, "[Log] IBFs construction finished\n");
    /* mark common sketches between alice and bob */
    for(i = 0; !err && i < alice.n; ++i) {
        if (!duplicate_alice[i]) {
            ibfbuf = (uint8_t*)1;/* Pure evil, but ok for a boolean variable */
            for(j = 0; !err && ibfbuf && j < bob.n; ++j) {
                if (!duplicate_bob[j]) {
                    a = &alice.a[i];
                    b = &bob.a[j];
                    if (!err) if (a->hash_width != b->hash_width || a->size != b->size || a->seed != b->seed) err = ERR_INCOMPATIBLE;
                    if (memcmp(alice.a[i].hashes, bob.a[i].hashes, bob.a[i].hash_width * bob.a[i].size) == 0) {
                        duplicate_alice[i] = 2;
                        duplicate_bob[j] = 2;
                        ibfbuf = NULL; /* using it as a bool variable, duplicate elements inside both alice and bob are already marked, so it is safe to end the inner loop */
                    }
                }
            }
        }
    }
    // fprintf(stderr, "[Log] exact cross matching done\n");
    /* ibf diff */
    diff.seres = NULL;
    diff.data = NULL;
    if (!err) err = ibf_sketch_diff(&ibf_alice, &ibf_bob, &diff);
    if (!err) err = ibf_sketch_store(opts.opath, &diff);
    callback_io.alice = &alice;
    callback_io.bob   = &bob;
    callback_io.dup_alice = duplicate_alice;
    callback_io.dup_bob   = duplicate_bob;
    if (!err) err = ibf_list_seq(&diff, get_sketches, &callback_io);
    // fprintf(stderr, "[Log] difference peeling done\n");
    /* check if all differences have been marked */
    if (!err) {
        j = 1;
        for(i = 0; i < alice.n; ++i) if (j && !duplicate_alice[i]) j = 0;
        if (!j) fprintf(stderr, "[Warning] Not all differences in Alice have been found\n");
        j = 1;
        for(i = 0; i < bob.n; ++i) if (j && !duplicate_bob[i]) j = 0;
        if (!j) fprintf(stderr, "[Warning] Not all differences in Bob have been found\n");
        // for(i = 0; i < alice.n; ++i) fprintf(stderr, "%u ", duplicate_alice[i]);
        // fprintf(stderr, "\n----------------------------\n");
        // for(i = 0; i < bob.n; ++i) fprintf(stderr, "%u ", duplicate_bob[i]);
        // fprintf(stderr, "\n");
    }
    assert(buflen + sizeof(diff.data->counter) == sizeof(bucket_t));
    if (!err) {
        fprintf(stdout, "%llu", diff.repetitions * (diff.chunk_size * (buflen + sizeof(diff.data->counter)) + sizeof(diff.seres->seed)) + 
                                sizeof(diff.repetitions) +
                                sizeof(diff.epsilon) +
                                sizeof(diff.chunk_size) 
                                #ifdef GLEN
                                + sizeof(diff.key_len)
                                #endif
                                );
    } else fprintf(stdout, "NaN");

    /* cleaning */
    if (!err) err = ibf_sketch_destroy(&diff);
    if (!err) err = ibf_sketch_destroy(&ibf_alice);
    if (!err) err = ibf_sketch_destroy(&ibf_bob);
    if (duplicate_alice) free(duplicate_alice);
    if (duplicate_bob) free(duplicate_bob);
    for(i = 0; !err && i < alice.n; ++i) if (!err) err = minhash_sketch_destroy(&alice.a[i]);
    kv_destroy(alice);
    for(i = 0; !err && i < bob.n; ++i) if (!err) err = minhash_sketch_destroy(&bob.a[i]);
    kv_destroy(bob);

    return err;
}

enum Error fill_minhash_sketch(param_t const * const opts, path_itr start, path_itr end, minhash_v_t * const mhdb) {
    gzFile fp;
    minhash_t dummy;
    kseq_t *seq;
    uint64_t i, j;
    struct stat st;
    enum Error err;
    assert(mhdb);
    fp = NULL;
    seq = NULL;
    err = NO_ERROR;
    dummy.size = 0;
    for(; !err && start != end; ++start) {
        kv_push(minhash_t, NULL, *mhdb, dummy);
        err = minhash_sketch_init(opts->s, opts->z, opts->w, &mhdb->a[mhdb->n-1]);
        if (strcmp(*start, "-") == 0) {
            if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
                fprintf(stderr, "Unable to use stdin as input\n");
                return ERR_FILE;
            }
        } else if (stat(*start, &st) == 0 && (st.st_mode & S_IFREG)) {
            if ((fp = gzopen(*start, "r")) == NULL) {
                fprintf(stderr, "Unable to open file %s\n", *start);
                return ERR_FILE;
            }
        } else return ERR_FILE;
        if (fp) seq = kseq_init(fp);
        while(!err && kseq_read(seq) >= 0) {
            for(i = 0, j = 0; !err && i < seq->seq.l; ++i) {
                if (seq_nt4_table[(uint8_t)seq->seq.s[i]] < 4) ++j;
                else j = 0;
                if (j >= opts->k) {
                    // minhash_print(stderr, &sketch_vec.a[sketch_vec.n-1]);
                    // fprintf(stderr, "i = %llu, j = %llu, kmer = %.*s\n", i, j, opts.k, &seq->seq.s[i-opts.k+1]);
                    err = minhash_sketch_insert(&seq->seq.s[i-opts->k+1], opts->k, opts->z, &mhdb->a[mhdb->n-1]);
                }
            }
        }
        if (seq) kseq_destroy(seq);
        if (fp) gzclose(fp);
        fp = NULL;
    }
    return err;
}

enum Error mark_duplicates(minhash_v_t const * const mhdb, uint8_t ** const duplicate) {
    uint64_t i, j;
    enum Error err;
    assert(mhdb);
    assert(*duplicate == NULL);
    err = NO_ERROR;
    if ((*duplicate = calloc(mhdb->n, sizeof **duplicate)) == 0) {
        fprintf(stderr, "Unable to allocate duplicate vector\n");
        err = ERR_ALLOC;
    }
    if (!err) {
        for(i = 0; i < mhdb->n; ++i) {
            for(j = i + 1; j < mhdb->n; ++j) {
                if (!(*duplicate)[j] && memcmp(mhdb->a[i].hashes, mhdb->a[j].hashes, mhdb->a[i].hash_width * mhdb->a[i].size) == 0) 
                    (*duplicate)[j] = 1;
            }
        }
    }
    return err;
}

void get_sketches(bucket_t const * const bucket, char origin, void *callback_io) {
    uint64_t i;
    callback_t *ino;
    assert(bucket);
    assert(callback_io);
    ino = (callback_t*)callback_io;
    /* mark differences */
    if (origin == 'i') {
        // fprintf(stderr, "found something unique to alice\n");
        for(i = 0; i < ino->alice->n; ++i) {
            if (!ino->dup_alice[i] && memcmp(ino->alice->a[i].hashes, bucket->keysum, ino->alice->a[i].hash_width * ino->alice->a[i].size) == 0)
                ino->dup_alice[i] = 3;
        }
    } else if (origin == 'j') {
        // fprintf(stderr, "found something unique to bob\n");
        for(i = 0; i < ino->bob->n; ++i) {
            if (!ino->dup_bob[i] && memcmp(ino->bob->a[i].hashes, bucket->keysum, ino->bob->a[i].hash_width * ino->bob->a[i].size) == 0)
                ino->dup_bob[i] = 3;
        }
    } else {
        fprintf(stderr, "[Warning] Unrecognized source ibf\n");
    }
    //return;
}

void print_collection_help() {
    fprintf(stderr, "[collection] options:\n");
    fprintf(stderr, "\t-o\tInvertible Bloom Filter file (binary output)\n");
    fprintf(stderr, "\t-a\tlist of fastx files making up Alice's database\n");
    fprintf(stderr, "\t-b\tlist of fastx files making up Bob's database\n");
    fprintf(stderr, "\t-k\tk-mer length [21]\n");
    fprintf(stderr, "\t-z\tnumber of hashes in each minHash sketch [1000]\n");
    fprintf(stderr, "\t-w\tbit-width of each hash [64]\n");
    fprintf(stderr, "\t-n\tupper bound on the number of different sketches to track (0 < n)\n");
    fprintf(stderr, "\t-r\tnumber of hash functions [3] (3 <= r <= 7)\n");
    fprintf(stderr, "\t-e\tepsilon [0] (0 <= epsilon <= 1)\n");
    fprintf(stderr, "\t-s\trandom seed [42]\n");
    fprintf(stderr, "\t-h\tshow this help\n");
    fprintf(stderr, "\nExample:\n");
    fprintf(stderr, "\tminhash_test -o <output_file> -a <fastx|gz> (*[fastx|gz]) -b <fastx|gz> (*[fastx|gz]) -k <k> -z <number of hashes> -w <hash width in bits> -s <seed>\n");
}

void init_options(param_t *parameters) {
    parameters->k = 21;
    parameters->z = 1000;
    parameters->s = 42;
    parameters->w = 64;
    parameters->n = 0;
    parameters->r = 3;
    parameters->e = 0;
    parameters->opath = NULL;
}

static enum Error option_parser(int *opt_idx, char opt, char** arg_ptr, param_t *parameters) {
    char *arg;
    static long int parsed_ld;
    static long long parsed_ll;
    if (arg_ptr != NULL) arg = *arg_ptr;
    else arg = NULL;
    switch(opt) {
        case 'o':
            parameters->opath = arg;
            break;
        case 'a':
            for(parameters->alice_start = parameters->alice_stop = arg_ptr; 
                *(parameters->alice_stop) != NULL && **(parameters->alice_stop) != '-'; 
                ++parameters->alice_stop) {++*opt_idx;}
            --*opt_idx;
            break;
        case 'b':
            for(parameters->bob_start = parameters->bob_stop = arg_ptr; 
                *(parameters->bob_stop) != NULL && **(parameters->bob_stop) != '-'; 
                ++parameters->bob_stop) {++*opt_idx;}
            --*opt_idx;
            break;
        case 'k':
            parsed_ld = strtol(arg, NULL, 10);
            if (parsed_ld > UINT8_MAX || parsed_ld < 0) return ERR_OUTOFBOUNDS;
            parameters->k = (unsigned char)parsed_ld;
            break;
        case 'z':
            parsed_ll = strtoll(arg, NULL, 10);
            if (parsed_ll > INT64_MAX || parsed_ll < 0) return ERR_OUTOFBOUNDS;
            parameters->z = (uint64_t)parsed_ll;
            break;
        case 'w':
            parsed_ld = strtol(arg, NULL, 10);
            if (parsed_ld > UINT8_MAX || parsed_ld < 0) return ERR_OUTOFBOUNDS;
            parameters->w = (unsigned char)parsed_ld;
            break;
        case 'n':
            parsed_ll = strtoll(arg, NULL, 10);
            if (parsed_ll > INT64_MAX || parsed_ll < 0) return ERR_OUTOFBOUNDS;
            parameters->n = (uint64_t)parsed_ll;
            break;
        case 'r':
            parsed_ld = strtol(arg, NULL, 10);
            if (parsed_ld > UINT8_MAX || parsed_ld < 0) return ERR_OUTOFBOUNDS;
            parameters->r = (unsigned char)parsed_ld;
            break;
        case 'e':
            parameters->e = atof(arg);
            break;
        case 's':
            parsed_ll = strtoll(arg, NULL, 10);
            if (parsed_ll > INT32_MAX || parsed_ll < 0) return ERR_OUTOFBOUNDS;
            parameters->s = (uint32_t)parsed_ll;
            break;
        case 'h':
            print_collection_help();
            break;
        default:
            return ERR_OPTION;
    }
    return NO_ERROR;
}

enum Error check_options(param_t const * const parameters) {
    if (parameters->k > 32) return ERR_OUTOFBOUNDS;
    if (parameters->w > 64) return ERR_OUTOFBOUNDS;
    if (parameters->n == 0) {
        fprintf(stderr, "Unspecified n"); 
        return ERR_OPTION;
    }
    if (parameters->e > 1 || parameters->e < 0) return ERR_OUTOFBOUNDS;
    if (parameters->r > 7 || parameters->r < 3) return ERR_OUTOFBOUNDS;
    return NO_ERROR;
}

/* Function for parsing arguments in the form -x :
 * 
 * @param argc Number of command line parameters
 * @param argv Vector of command line parameters in string format
 * @opt_string argparse-like string listing the available options with : indicating if argument is required
 * @parser Function for parsing each option in opt_string, saving the parsed values inside params
 * @params Struct containing the parsed arguments.
 * @return 0 - success, 1 - unknown option, any error code returned by the parser function
 */
enum Error parse_options(int argc, char *argv[], char const * const opt_string, enum Error (*parser)(int*, char, char**, param_t*), param_t *params) {
    int i, j, ok;
    //int opt_ind;
    char **opt_arg;
    char short_opt;
    size_t opt_len;

    assert(argv);
    assert(opt_string);
    assert(parser);
    assert(params);

    //opt_ind = -1;
    opt_arg = NULL;
    short_opt = '\0';
    for(i = 1, ok = 1; i < argc && ok; ++i) {
        // fprintf(stderr, "now at %s\n", argv[i]);
        if (argv[i][0] == '-') {/*Order is important!*/
            opt_len = strlen(argv[i]);
            if(opt_len == 1) short_opt = argv[i][0];/* just '-' if everything fine */
            else if (opt_len == 2) short_opt = argv[i][1];
            else ok = 0;/* long option not supported */
            opt_len = strlen(opt_string);
            for(j = 0; j < opt_len; ++j)
                if (short_opt == opt_string[j] && short_opt != ':') break;

            if (j < opt_len) { /* short_opt was recognized */
                if (opt_string[j + 1] == ':') { /* if option requires argument then set variables */
                    //opt_ind = 
                    ++i;
                    opt_arg = &argv[i];
                }
                j = parser(&i, short_opt, opt_arg, params);
                if (j) return j;
            } else { /* unknown option */
                return ERR_OPTION;
            }
        } else {
            ok = 0;
            --i;
        }
    }
    params->opt_end = i;
    return NO_ERROR;
}

/*
    if ((duplicate_alice = calloc(alice.n, sizeof *duplicate_alice)) == 0) {
        fprintf(stderr, "Unable to allocate duplicate vector\n");
        err = ERR_ALLOC;
    }
    if ((duplicate_bob = calloc(bob.n, sizeof *duplicate_bob)) == 0) {
        fprintf(stderr, "Unable to allocate duplicate vector\n");
        err = ERR_ALLOC;
    }
    if (!err) {
        for(i = 0; i < alice.n; ++i) {
            for(j = i + 1; j < alice.n; ++j) {
                if (!duplicate_alice[j] && memcmp(alice.a[i].hashes, alice.a[j].hashes, alice.a[i].hash_width * alice.a[i].size) == 0) 
                    duplicate_alice[j] = 1;
            }
        }
    }
    if (!err) {
        for(i = 0; i < bob.n; ++i) {
            for(j = i + 1; j < bob.n; ++j) {
                if (!duplicate_bob[j] && memcmp(bob.a[i].hashes, bob.a[j].hashes, bob.a[i].hash_width * bob.a[i].size) == 0) 
                    duplicate_bob[j] = 1;
            }
        }
    }
*/
