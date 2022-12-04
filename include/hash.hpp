#ifndef HASH_HPP
#define HASH_HPP

#include "MurmurHash3.hpp"

namespace hash {

template <typename T>
static inline uint64_t hash64(T val, uint64_t seed)
{
    uint64_t hval[2];
    MurmurHash3_x64_128(reinterpret_cast<void*>(&val), sizeof(T), seed, hval);
    return hval[0];
}

static inline uint64_t hash64(uint8_t const* key, uint32_t keylen, uint32_t seed)
{
    uint64_t hash;
    assert(keylen < std::numeric_limits<int>::max());
    MurmurHash3_x64_128(reinterpret_cast<const void*>(key), keylen, seed, reinterpret_cast<void*>(&hash));
    return hash;
}

class mm_pos_extractor {
    public:

        std::size_t operator()(uint64_t kmer) {
            uint64_t mval = hash64<uint64_t>(kmer & mask, 0);
            uint8_t minpos = 0;
            for (std::size_t i = 0; i < k - m + 1; ++i) {
                auto val = hash64<uint64_t>(kmer & mask, 0);
                if (mval >= val) {
                    mval = val;
                    minpos = i;
                }
                kmer >>= 2;
            }
            return k - m - minpos; // CHECK
        }
    private:
        uint8_t k;
        uint8_t m;
        uint64_t mask;
};

}

#endif // HASH_HPP