#include "endian_fixer.h"

static const int num = 42;/* The answer is 42 */

uint64_t hton64(uint64_t value) {
    if (*(const char*)&num == num) { /* little endian */
        const uint32_t high_part = htonl((uint32_t)(value >> 32));
        const uint32_t low_part = htonl((uint32_t)(value & 0xFFFFFFFFLL));
        return ( ((uint64_t)low_part) << 32 ) | high_part;
    } else { /* big endian */
        return value;
    }
}
