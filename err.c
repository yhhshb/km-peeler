#include <stdio.h>
#include "err.h"

enum Error translate_error(enum Error code) {
    switch (code) {
        case ERR_OPTION:
            fprintf(stderr, "[OptionError]");
            break;
        case ERR_FILE:
            fprintf(stderr, "[FileError]");
            break;
        case ERR_IO:
            fprintf(stderr, "[IOError]");
            break;
        case ERR_VALUE:
            fprintf(stderr, "[ValueError]");
            break;
        case ERR_ALLOC:
            fprintf(stderr, "[AllocationError]");
            break;
        case ERR_OUTOFBOUNDS:
            fprintf(stderr, "[OutOfBoundsError]");
            break;
        case ERR_RUNTIME:
            fprintf(stderr, "[RuntimeError]");
            break;
        case ERR_INCOMPATIBLE:
            fprintf(stderr, "[IncompatibleError]");
            break;
        default:
            fprintf(stderr, "[UnrecognisedError]");
    }
    return code;
}

enum Error print_error(enum Error code, char *msg) {
    if (code != NO_ERROR) {
        translate_error(code);
        fprintf(stderr, " %s\n", msg);
    }
    return code;
}
