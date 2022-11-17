#include "utils.hpp"
#include "MurmurHash3.hpp"
#include <cassert>
#include <limits>

uint64_t hash64(uint8_t const* key, uint32_t keylen, uint32_t seed)
{
    uint64_t hash;
    assert(keylen < std::numeric_limits<int>::max());
    MurmurHash3_x64_128(reinterpret_cast<const void*>(key), keylen, seed, reinterpret_cast<void*>(&hash));
    return hash;
}

void shift_byte_vec(uint8_t* dst, uint8_t* src, std::size_t size, uint8_t shift)
{
    assert(shift < 64);
    if (not size) return;
    for (std::size_t i = 0; i < size-1; ++i,++src,++dst)
    {
        *dst = ((*src)<<shift) | (*(src+1)>>(8-shift));
    }
    *dst = ((*src)<<shift);
}