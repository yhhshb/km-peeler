#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "kmers_main.h"
#include "minimizers_main.h"
#include "syncmers_main.h"
#include "sample_main.h"
#include "build_main.h"
#include "diff_main.h"
#include "list_main.h"
#include "jaccard_main.h"
#include "collection_main.h"
#include "print_main.h"
#include "dump_main.h"
#include "ketopt.h"
#include "err.h"

void print_subcommands();

int main(int argc, char** argv) {
    ketopt_t om = KETOPT_INIT;
    int c;
    enum Error error_code;
    while((c = ketopt(&om, argc, argv, 0, "h", 0)) >= 0)
    {
        if(c == 'h')
        {
            print_subcommands();
            return EXIT_SUCCESS;
        }
    }
    if(om.ind == argc) {
        fprintf(stderr, "[Error] Expecting subcommand: ");
        for(c = 0; c < argc; ++c) fprintf(stderr, "%s ", argv[c]);
        fprintf(stderr, "\n");
        return EXIT_FAILURE;
    }
    if (strcmp(argv[om.ind], "kmers") == 0) {
        error_code = kmers_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "minimizers") == 0) {
        error_code = minimizers_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "syncmers") == 0) {
        error_code = syncmers_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "sample") == 0) {
        error_code = sample_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "build") == 0) {
        error_code = build_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "diff") == 0) {
        error_code = diff_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "list") == 0) {
        error_code = list_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "jaccard") == 0) {
        error_code = jaccard_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "collection") == 0) {
        error_code = collection_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "count") == 0) {
        error_code = count_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "print") == 0) {
        error_code = print_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "dump") == 0) {
        error_code = dump_main(argc - om.ind, &argv[om.ind]);
    } else if (strcmp(argv[om.ind], "WSIZE") == 0) {
        fprintf(stderr, "%llu B\n", WSIZE);
        error_code = NO_ERROR;
    } else {
        fprintf(stderr, "Missing subcommand\n\n");
        print_subcommands();
        return EXIT_FAILURE;
    }

    if (error_code != NO_ERROR) {
        fprintf(stderr, "[FAILURE] Error: ");
        print_error(error_code, "");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void print_subcommands() {
    fprintf(stderr, "ibfseq [subcommand] [options]\n");
    fprintf(stderr, "\tavailable subcommands are:\n");
    fprintf(stderr, "\tkmers\tfragment sequences into k-mers\n");
    fprintf(stderr, "\tminimizers\tfragment sequences by using minimizers\n");
    fprintf(stderr, "\tsyncmers\tfragment sequences by using syncmers\n");
    fprintf(stderr, "\tsample\tsample sequences\n");
    fprintf(stderr, "\tbuild\tInvertible Bloom Filter construction\n");
    fprintf(stderr, "\tdiff\tcompute difference between two Invertible Bloom Filters\n");
    fprintf(stderr, "\tlist\ttry to list the content of an Invertible Bloom Filter\n");
    fprintf(stderr, "\tjaccard\tcompute jaccard similarity between two Invertible Bloom Filters\n");
    fprintf(stderr, "\tcollection\tbuild an IBF storing a collection of minHash sketches\n");
    fprintf(stderr, "\tcount\tcount elements inside IBF\n");
    fprintf(stderr, "\tprint\tprint the bucket counts of an Invertible Bloom Filter\n");
    fprintf(stderr, "Use option -h of a subcommand to show its options\n");
}
