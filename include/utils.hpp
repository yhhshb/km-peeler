#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <limits>

namespace kmp {

template <typename IntegerType1, typename IntegerType2>
IntegerType1 ceil_size(IntegerType1 x, IntegerType2 y)
{
    return (x + y - 1) / y;
}

void shift_byte_vec(uint8_t* dst, uint8_t* src, std::size_t size, uint8_t shift);

} // namespace kmp

#endif // UTILS_HPP