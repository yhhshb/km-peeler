#include <string>
#include <vector>
#include <cassert>

#include "../include/utils.hpp"
#include "../include/constants.hpp"

namespace kmp {

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

std::string vec2kmer(std::vector<uint8_t> const& packed, uint8_t k)
{
    assert(ceil_size(2*k, 8) == packed.size());
    std::string ret(k, 'N');
    for (std::size_t i = 0; i < packed.size(); ++i) {
        uint8_t b = packed.at(i);
        int lb = i == 0 ? (packed.size() * 4 - k) : 4;
        assert(lb >= 0);
        for (int j = 4; j > lb; --j) {
            ret[i + j] = packed[i] & 0x03;
            b >>= 2;
        }
    }
}

} // namespace kmp