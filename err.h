#ifndef ERR_H
#define ERR_H

#define NO_ERROR 0
enum Error {ERR_OPTION = 1, ERR_FILE, ERR_IO, ERR_VALUE, ERR_ALLOC, ERR_OUTOFBOUNDS, ERR_RUNTIME, ERR_INCOMPATIBLE};

enum Error translate_error(enum Error code);

enum Error print_error(enum Error code, char *msg);

#endif/*ERR_H*/
