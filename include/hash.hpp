#ifndef HASH_HPP
#define HASH_HPP

#include "MurmurHash3.hpp"

namespace hash {

class double_hash64
{
    public:
        typedef uint64_t hash_type;

        std::array<uint64_t, 2> operator()(uint8_t const* key, uint32_t keylen, uint32_t seed) const noexcept
        {
            std::array<uint64_t, 2> hval;
            assert(keylen < std::numeric_limits<int>::max());
            MurmurHash3_x64_128(reinterpret_cast<const void*>(key), keylen, seed, hval.data());
            return hval;
        }

        template <typename T>
        std::array<uint64_t, 2> operator()(T val, uint64_t seed) const noexcept
        {
            return operator()(reinterpret_cast<uint8_t*>(&val), sizeof(T), seed);
        }
};

class hash64
{
    public:
        typedef uint64_t hash_type;

        uint64_t operator()(uint8_t const* key, uint32_t keylen, uint32_t seed) const noexcept
        {
            return hasher(key, keylen, seed)[0];
        }

        template <typename T>
        uint64_t operator()(T val, uint64_t seed) const noexcept
        {
            return operator()(reinterpret_cast<uint8_t*>(&val), sizeof(T), seed);
        }

    private:
        double_hash64 hasher;
};

class mm_pos_extractor {
    public:
        std::size_t operator()(uint64_t kmer) {
            uint64_t mval = hasher(kmer & mask, 0);
            uint8_t minpos = 0;
            for (std::size_t i = 0; i < k - m + 1; ++i) {
                auto val = hasher(kmer & mask, 0);
                if (mval >= val) {
                    mval = val;
                    minpos = i;
                }
                kmer >>= 2;
            }
            return k - m - minpos; // TODO check correctness
        }
    private:
        uint8_t k;
        uint8_t m;
        uint64_t mask;
        hash64 hasher;
};

}

#endif // HASH_HPP