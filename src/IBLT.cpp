#include "IBLT.hpp"
#include "utils.hpp"
#include <cmath>

namespace kmp {

IBLT::IBLT(uint8_t k, uint8_t r, double epsilon, uint64_t n, uint64_t seed = 42) : klen(k), reps(r), eps(epsilon), max_diff(n), pseed(seed)
{
    if (reps < 3 or reps > 7) throw std::invalid_argument("peelability constants are defined only for 3 <= r <= 7");
    uint64_t chunk_size = (ck_table[reps] + eps) * max_diff / reps + 1;
    number_of_buckets = chunk_size * reps;
    hf_size = static_cast<decltype(hf_size)>(std::ceil((reps - 2) * std::log(static_cast<double>(max_diff)) + reps));
    bucket_size = kmp::ceil(2 * uint64_t(klen) + uint64_t(hf_size), 8ULL); // hash and payload packed together;
    counts.resize(kmp::ceil(number_of_buckets, 4ULL)); // 4 counts per byte
    hp_buckets.resize(number_of_buckets * bucket_size);
}

void IBLT::insert(std::vector<uint8_t> const& kmer) noexcept
{
    assert(kmer.size() >= klen / 4);
    for (uint8_t i = 0; i < reps; ++i) {
        std::size_t idx = hash64(kmer.data(), kmer.size(), i) % (number_of_buckets / reps);
        inc_count_at(idx);
        xor_at(idx, kmer);
    }
}

void IBLT::remove(std::vector<uint8_t> const& kmer) noexcept
{
    assert(kmer.size() >= klen / 4);
    for (uint8_t i = 0; i < reps; ++i) {
        std::size_t idx = hash64(kmer.data(), kmer.size(), i) % (number_of_buckets / reps);
        dec_count_at(idx);
        xor_at(idx, kmer);
    }
}

std::vector<std::vector<uint8_t>> IBLT::list()
{
    std::vector<std::vector<uint8_t>> res;
    std::size_t idx;
    while ((idx = find_peelable_bucket()) != std::numeric_limits<decltype(idx)>::max()) {
        res.push_back(get_payload(idx));
        remove(res.back());
    }
    return res;
}

void IBLT::update_count_at(std::size_t idx, std::function<uint8_t(uint8_t)> op) noexcept
{
    std::size_t byte_idx = idx / 4;
    std::size_t bit_couple_idx = idx % 4;
    uint8_t shift = (3 - bit_couple_idx) * 2;
    uint8_t c = (counts.at(byte_idx) >> shift) & 0x03;
    counts[byte_idx] &= ~(0x03 << shift);
    c = op(c);
    c &= 0x03;
    counts[byte_idx] |= c << shift;
}

void IBLT::inc_count_at(std::size_t idx) noexcept
{
    update_count_at(idx, [](uint8_t c){return c + 1;});
}

void IBLT::dec_count_at(std::size_t idx) noexcept
{
    update_count_at(idx, [](uint8_t c){return c - 1;});
}

void IBLT::xor_at(std::size_t idx, std::vector<uint8_t> const& kmer) noexcept
{
    
}

std::size_t IBLT::find_peelable_bucket() const
{

}

std::vector<uint8_t> IBLT::get_payload(std::size_t idx) const noexcept
{

}

std::ostream& operator<<(std::ostream& os, const IBLT& obj)
{
    // write obj to stream
    return os;
}

}
