#ifndef UTILS_HPP
#define UTILS_HPP

namespace kmp {

template <typename IntegerType>
IntegerType ceil(IntegerType x, IntegerType y)
{
    return (x + y - 1) / y;
}

uint64_t hash64(uint8_t const* key, uint64_t keylen, uint32_t seed);

void shift_byte_vec(uint8_t* dst, uint8_t* src, std::size_t size, uint8_t shift);

} // namespace kmp

#endif // UTILS_HPP