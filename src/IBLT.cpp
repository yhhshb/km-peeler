#include "IBLT.hpp"
#include "hash.hpp"
#include "utils.hpp"
#include <cmath>

namespace kmp {

#ifndef NDEBUG
const std::array<uint8_t, 8> left_align_mask = {0x00, 0x80, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC, 0xFE};
#endif

IBLT::IBLT(uint8_t k, uint8_t r, double epsilon, uint64_t n, uint64_t seed = 42) 
    : 
    klen(k), 
    reps(r), 
    eps(epsilon), 
    max_diff(n), 
    pseed(seed),
    hrc_bit_size(static_cast<decltype(hrc_bit_size)>(std::ceil((reps - 2) * std::log(static_cast<double>(max_diff)) + reps))),
    prefix_len(hrc_bit_size % 8),
    mask(~((uint8_t(1) << (8-prefix_len)) - 1)),
    number_of_inserted_items(0)
{
    assert(mask == left_align_mask[prefix_len]);
    if (reps < 3 or reps > 7) throw std::invalid_argument("peelability constants are defined only for 3 <= r <= 7");
    uint64_t chunk_size = (ck_table[reps] + eps) * max_diff / reps + 1;
    number_of_buckets = chunk_size * reps;
    if (hrc_bit_size > sizeof(decltype(hasher)::hash_type) * 8) throw std::runtime_error("[IBLT] Bucket check hash size larger than hash function return type");
    bucket_size = kmp::ceil_size(2ULL * klen + hrc_bit_size, 8); // hash and payload packed together;
    counts.resize(kmp::ceil_size(number_of_buckets, 4)); // 4 count values in each byte
    hp_buckets.resize(number_of_buckets * bucket_size);
    std::size_t hrc_byte_size = kmp::ceil_size(hrc_bit_size, 8);
    hash_redundancy_code_buffer.resize(hrc_byte_size);
    payload_buffer.resize(kmp::ceil_size(2ULL * klen, 8));
}

void IBLT::insert(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept
{
    assert(kmer_byte_size >= klen / 4);
    for (uint8_t i = 0; i < reps; ++i) {
        auto hval = hasher(kmer, kmer_byte_size, i);
        std::size_t idx = hval[0] % (number_of_buckets / reps) + (number_of_buckets / reps) * i;
        inc_count_at(idx);
        xor_at(idx, kmer, hval[1]);
    }
    ++number_of_inserted_items;
}

void IBLT::remove(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept
{
    assert(kmer_byte_size >= klen / 4);
    for (uint8_t i = 0; i < reps; ++i) {
        auto hval = hasher(kmer, kmer_byte_size, i);
        std::size_t idx = hval[0] % (number_of_buckets / reps) + (number_of_buckets / reps) * i;
        dec_count_at(idx);
        xor_at(idx, kmer, hval[1]);
    }
    --number_of_inserted_items;
}

void IBLT::subtract(IBLT const& other)
{
    constexpr bool hashers_equality = std::is_same<decltype(hasher), decltype(other.hasher)>::value;
    if (
        klen != other.klen or
        reps != other.reps or
        max_diff != other.max_diff or
        pseed != other.pseed or
        number_of_buckets != other.number_of_buckets or
        not hashers_equality or
        counts.size() != other.counts.size() or
        hp_buckets.size() != other.hp_buckets.size()
    ) throw std::runtime_error("[IBLT::subtract] Incompatible IBLTs");

    for (std::size_t i = 0; i < hp_buckets.size(); ++i) hp_buckets[i] ^= other.hp_buckets[i]; // FIXME optimize
    // TODO how to efficiently subtract count remainders ?

    std::fill(std::begin(payload_buffer), std::end(payload_buffer), 0);
    std::fill(std::begin(hash_redundancy_code_buffer), std::end(hash_redundancy_code_buffer), 0);

    if (number_of_inserted_items >= other.number_of_inserted_items) number_of_inserted_items -= other.number_of_inserted_items;
    else number_of_inserted_items = other.number_of_inserted_items - number_of_inserted_items;
}

IBLT::failure_t IBLT::list(std::vector<std::vector<uint8_t>>& res) noexcept
{
    idx_stack_t peelable_idx_stack;
    std::size_t tmp;
    while ((tmp = find_peelable_bucket()) != number_of_buckets && res.size() < number_of_inserted_items + 1) {
        peelable_idx_stack.push(tmp);
        while (not peelable_idx_stack.empty()) {
            auto idx = peelable_idx_stack.top();
            peelable_idx_stack.pop();
            res.push_back(get_payload(idx));
            peel(res.back().data(), res.size(), idx, peelable_idx_stack);
        }
    }
    if (res.size() > number_of_inserted_items) return INFINITE_LOOP;
    if (res.size() < number_of_inserted_items) return UNPEELABLE;
    return NONE;
}

std::size_t IBLT::size() const noexcept
{
    return number_of_inserted_items;
}

void IBLT::inc_count_at(std::size_t idx) noexcept
{
    update_count_at(idx, [](uint8_t c){return c + 1;});
}

void IBLT::dec_count_at(std::size_t idx) noexcept
{
    update_count_at(idx, [](uint8_t c){return c - 1;});
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

uint8_t IBLT::get_count_at(std::size_t idx) const noexcept
{
    std::size_t byte_idx = idx / 4;
    std::size_t bit_couple_idx = idx % 4;
    uint8_t shift = (3 - bit_couple_idx) * 2;
    return (counts.at(byte_idx) >> shift) & 0x03;
}

void IBLT::xor_at(std::size_t idx, uint8_t const * const kmer, uint64_t hrc) noexcept // hrc = hash redundancy code
{
    const std::size_t bucket_byte_start = idx * bucket_size;
    std::size_t payload_shift = hrc_bit_size / 8;

    assert(kmer);
    assert(payload_shift <= sizeof(decltype(hrc)));

    uint8_t *const hrc_ptr = reinterpret_cast<uint8_t*>(&hrc);
    // XOR hash value of k-mer to bucket (left justified) except the last byte (if any)
    for (std::size_t i = 0; i < payload_shift; ++i) hp_buckets[i + bucket_byte_start] ^= hrc_ptr[i]; 
    if (prefix_len) { // if the last byte is packed, mask and XOR
        
        hrc_ptr[payload_shift] &= mask;
        hp_buckets[payload_shift + bucket_byte_start] ^= hrc_ptr[payload_shift];
        // if (klen != 0 && (mask & kmer[0]) != 0) throw std::runtime_error("[IBLT::xor_at] Hash check and k-mer overlap one another");
        assert(klen == 0 or (mask & kmer[0]) == 0);
    } else { // here everything is byte-aligned
        ++payload_shift;
    }
    // XOR k-mer to bucket (right-justified)
    for (std::size_t i = 0; i < klen/4; ++i) hp_buckets[bucket_byte_start + payload_shift + i] ^= kmer[i]; // klen/4 ok since payload_shift takes care of packed byte
    assert(payload_shift + klen/4 == bucket_size);
}

std::size_t IBLT::unpack_at(std::size_t idx) noexcept
{
    const std::size_t bucket_byte_start = idx * bucket_size;
    std::size_t payload_shift = hrc_bit_size / 8;

    auto start_itr = hp_buckets.cbegin() + bucket_byte_start;
    hash_redundancy_code_buffer.assign(start_itr, start_itr + hash_redundancy_code_buffer.size());
    if (prefix_len) hash_redundancy_code_buffer.back() &= mask; // next: start copying payload from same end byte of hrc
    else ++payload_shift; // otherwise, if everything was perfectly aligned, start from next byte
    payload_buffer.assign(start_itr + payload_shift, start_itr + bucket_size);
    if (prefix_len) payload_buffer.front() &= ~mask;
    // TODO Find good error checking for this part
    assert(number_of_buckets % reps == 0);
    return idx / (number_of_buckets / reps);
}

bool IBLT::is_peelable(std::size_t idx) noexcept
{
    auto c = get_count_at(idx);
    if (c == 0 or c == 2) return false;
    std::size_t row = unpack_at(idx); // fill buffers with bucket values
    decltype(hasher)::hash_type cval = hasher(payload_buffer.data(), payload_buffer.size(), row)[1];
    auto cval_ptr = reinterpret_cast<uint8_t*>(&cval);
    cval_ptr[hash_redundancy_code_buffer.size()-1] &= mask;
    return std::equal(hash_redundancy_code_buffer.cbegin(), hash_redundancy_code_buffer.cend(), cval_ptr); // check hrc and recomputed hash value
}

std::size_t IBLT::find_peelable_bucket() noexcept
{
    for (std::size_t i = 0; i < number_of_buckets; ++i) {
        if (is_peelable(i)) return i;
    }
    return number_of_buckets;
}

std::vector<uint8_t> IBLT::get_payload(std::size_t idx) noexcept
{
    [[maybe_unused]] unpack_at(idx);
    return payload_buffer;
}

void IBLT::peel(uint8_t const * const kmer, std::size_t kmer_byte_size, std::size_t origin_bucket, idx_stack_t& idxs)
{ // this is basically a modified version of remove. New peelable buckets are added to the stack if they are found.
    assert(kmer_byte_size >= klen / 4);
    for (uint8_t i = 0; i < reps; ++i) {
        auto hval = hasher(kmer, kmer_byte_size, i);
        std::size_t idx = hval[0] % (number_of_buckets / reps) + (number_of_buckets / reps) * i;
        dec_count_at(idx);
        xor_at(idx, kmer, hval[1]);
        if (is_peelable(idx)) idxs.push(idx);
    }
}

std::ostream& operator<<(std::ostream& os, const IBLT& obj)
{
    // write obj to stream
    return os;
}

}
