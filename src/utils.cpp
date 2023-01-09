#include "../include/utils.hpp"
#include <cassert>

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