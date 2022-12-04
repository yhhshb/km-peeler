#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <limits>

namespace kmp {

template <typename IntegerType>
IntegerType ceil(IntegerType x, IntegerType y)
{
    return (x + y - 1) / y;
}

void shift_byte_vec(uint8_t* dst, uint8_t* src, std::size_t size, uint8_t shift);

} // namespace kmp

#endif // UTILS_HPP