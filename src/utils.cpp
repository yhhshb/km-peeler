#include <string>
#include <array>
#include <vector>
#include <sstream>
#include <cassert>

#include "../include/utils.hpp"
#include "../bundled/biolib/include/constants.hpp"

#include <iostream>

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
    assert(std::size_t(ceil_size(2*k, 8)) == packed.size());
    std::stringstream strm;
    for (std::size_t i = 0; i < packed.size(); ++i) {
        uint8_t b = packed.at(i);
        for (int j = 0; j < 4; ++j) {
            strm << constants::bases[(b & 0xC0) >> 6];
            b <<= 2;
        }
    }
    return strm.str().substr(packed.size()*4-k, k);
}

void print_byte_vec(std::ostream& os, uint8_t const * const v, std::size_t vlen)
{
    for (std::size_t i = 0; i < vlen; ++i) {
        for (std::size_t j = 7; j != std::numeric_limits<std::size_t>::max(); --j) {
            os << uint32_t((v[i] >> j) & 0x01);
        }
        os << " ";
    }
}

void little2big(uint8_t* val_bytes, std::size_t val_size)
{
    for (std::size_t i = 0; i < val_size / 2; i++)
    {
        uint8_t temp = val_bytes[i];
        val_bytes[i] = val_bytes[val_size - 1 - i];
        val_bytes[val_size - 1 - i] = temp;
    }
}

void dump_byte_vec(std::FILE* outstrm, uint8_t const * const val_bytes, std::size_t val_size) 
{
    for(std::size_t i = 0; i < val_size; ++i) {
		fprintf(outstrm, "%02X", val_bytes[i]);
	}
}

} // namespace kmp