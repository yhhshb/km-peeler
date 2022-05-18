#ifndef ENDIAN_FIXER_H
#define ENDIAN_FIXER_H

#include <stdint.h>
#include <arpa/inet.h>

#define hton8(v) (v)
#define ntoh8(v) (hton8(v))
#define hton16(v) (htons(v))
#define ntoh16(v) (ntohs(v))
#define hton32(v) (htonl(v))
#define ntoh32(v) (ntohl(v))

uint64_t hton64(uint64_t value);
#define ntoh64(v) (hton64(v))

#endif/*ENDIAN_FIXER_H*/
