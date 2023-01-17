#include "../include/IBLT.hpp"

#include <cmath>
#include "../include/io.hpp"
#include "../include/utils.hpp"

namespace kmp {

#ifndef NDEBUG
const std::array<uint8_t, 8> left_align_mask = {0x00, 0x80, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC, 0xFE};
#endif

IBLT::IBLT() : klen(0), reps(0), eps(0), max_diff(0), pseed(0), hrc_bit_size(0), prefix_len(0), mask(0), number_of_inserted_items(0)
{}

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
    resize_buffers();
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
    auto subtract_count_bytes = [](uint8_t a, uint8_t b) {
        uint8_t diff = 0;
        for (uint8_t i = 0; i < 4; ++i) {
            diff |= ((uint8_t(a & 0x03) - uint8_t(b & 0x03)) & 0x03) << (2*i);
            a >>= 2;
            b >>= 2;
        }
        return diff;
    };
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
    for (std::size_t i = 0; i < counts.size(); ++i) counts[i] = subtract_count_bytes(counts.at(i), other.counts.at(i)); // FIXME optimize

    std::fill(std::begin(payload_buffer), std::end(payload_buffer), 0); // reset buffers, just to be sure
    std::fill(std::begin(hash_redundancy_code_buffer), std::end(hash_redundancy_code_buffer), 0);
    if (number_of_inserted_items >= other.number_of_inserted_items) number_of_inserted_items -= other.number_of_inserted_items;
    else number_of_inserted_items = other.number_of_inserted_items - number_of_inserted_items;
}

IBLT::failure_t IBLT::list(std::vector<std::vector<uint8_t>>& positives, std::vector<std::vector<uint8_t>>& negatives)
{ // FIXME carefully check for errors during peeling
    idx_stack_t peelable_idx_stack;
    std::size_t tmp;
    const std::size_t peeling_bound = 2 * max_diff;
    while (
        (tmp = find_peelable_bucket()) != number_of_buckets and 
        (positives.size() + negatives.size()) < peeling_bound) // FIXME
    {
        peelable_idx_stack.push(tmp);
        while (not peelable_idx_stack.empty()) {
            auto idx = peelable_idx_stack.top();
            peelable_idx_stack.pop();
            unpack_at(idx);
            if (idx >= 0) positives.push_back(payload_buffer);
            else negatives.push_back(payload_buffer);
            peel(payload_buffer.data(), payload_buffer.size(), idx, peelable_idx_stack);
        }
    }
    std::size_t peeled_count = positives.size() + negatives.size();
    if (peeled_count >= peeling_bound) return INFINITE_LOOP;
    if (peeled_count < peeling_bound) return UNPEELABLE;
    return NONE;
}

IBLT::failure_t IBLT::list(std::size_t& positive_size, std::size_t& negative_size) 
{
    std::vector<std::vector<uint8_t>> p, n;
    auto ret = list(p, n);
    positive_size = p.size();
    negative_size = n.size();
    return ret;
}

std::size_t IBLT::size() const noexcept
{
    return number_of_inserted_items;
}

uint8_t IBLT::get_k() const noexcept
{
    return klen;
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

void IBLT::resize_buffers()
{
    std::size_t hrc_byte_size = kmp::ceil_size(hrc_bit_size, 8);
    hash_redundancy_code_buffer.resize(hrc_byte_size);
    payload_buffer.resize(kmp::ceil_size(2ULL * klen, 8));
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
    unpack_at(idx);
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

IBLT load(std::string filename, std::size_t& byte_size)
{
    IBLT iblt;
    byte_size = io::load(iblt, filename);
    iblt.resize_buffers();
    return iblt;
}

std::ostream& operator<<(std::ostream& os, const IBLT& obj)
{
    // write obj to stream
    return os;
}

std::ostream& operator<<(std::ostream& os, IBLT::failure_t const& err)
{
    switch(err) {
        case IBLT::UNPEELABLE:
            os << "None";
            break;
        case IBLT::INFINITE_LOOP:
            os << "inf";
            break;
        default:
            throw std::runtime_error("[Formatting] Unrecognized IBLT error code");
    }
}

}
