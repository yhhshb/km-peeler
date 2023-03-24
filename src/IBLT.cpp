#include <cmath>
#include <cassert>

#include "../include/IBLT.hpp"
#include "../include/io.hpp"
#include "../include/utils.hpp"

#ifndef NDEBUG
#include <iostream>
#endif

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
    hrc_bit_size(static_cast<decltype(hrc_bit_size)>(std::ceil((reps - 2) * std::log2(static_cast<double>(max_diff)) + reps))),
    prefix_len(hrc_bit_size % 8),
    mask(~((uint8_t(1) << (8-prefix_len)) - 1)),
    number_of_inserted_items(0)
{
    assert(mask == left_align_mask[prefix_len]);
    if (reps < 3 or reps > 7) throw std::invalid_argument("peelability constants are defined only for 3 <= r <= 7");
    uint64_t chunk_size = (ck_table[reps] + eps) * max_diff / reps + 1;
    number_of_buckets = chunk_size * reps;
    if (number_of_buckets > std::size_t(std::numeric_limits<long long>::max())) throw std::runtime_error("[IBLT] maximum number of buckets is 2^63 - 1 (signed long)");
    if (hrc_bit_size > sizeof(decltype(hasher)::hash_type) * 8) throw std::runtime_error("[IBLT] Bucket check hash size larger than hash function return type");
    bucket_size = kmp::ceil_size(2ULL * klen + hrc_bit_size, 8); // hash and payload packed together;
    counts.resize(kmp::ceil_size(number_of_buckets, 4)); // 4 count values in each byte
    hp_buckets.resize(number_of_buckets * bucket_size);
    resize_buffers();
}

void IBLT::insert(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept
{// TODO deal with endianness
    assert(payload_buffer.size());
    assert(kmer_byte_size >= payload_buffer.size());
    for (uint8_t i = 0; i < reps; ++i) {
        auto [idx, cval] = hash_kmer(kmer + kmer_byte_size - payload_buffer.size(), i);
        inc_count_at(idx);
        xor_at(idx, kmer + kmer_byte_size - payload_buffer.size(), cval);
    }
    ++number_of_inserted_items;
}

void IBLT::remove(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept
{
    assert(payload_buffer.size());
    assert(kmer_byte_size >= payload_buffer.size());
    for (uint8_t i = 0; i < reps; ++i) {
        auto [idx, cval] = hash_kmer(kmer + kmer_byte_size - payload_buffer.size(), i);
        dec_count_at(idx);
        xor_at(idx, kmer + kmer_byte_size - payload_buffer.size(), cval);
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
{
    auto find_peelable_bucket = [this]()
    {
        bool negative = false;
        for (std::size_t i = 0; i < number_of_buckets; ++i)
            if (is_peelable(i, negative)) 
                return peel_bucket_state_t(i, negative);
        return peel_bucket_state_t();
    };
    const std::size_t peeling_max_iterations = 2 * max_diff;
    std::size_t peeled_count = 0;
    std::vector<uint8_t> buffer;
    buffer.resize(payload_buffer.size());
    peel_bucket_state_t next;
    while ((next = find_peelable_bucket()).is_valid() and peeled_count < peeling_max_iterations)
    {
        while (next.is_valid() and peeled_count < peeling_max_iterations) {
            assert(buffer.size() == payload_buffer.size());
            buffer = payload_buffer; // make local copy of kmer since payload_buffer is a shared "global" method variable
            unpack_at(next.get_index());
            if (next.is_negative()) negatives.push_back(payload_buffer);
            else positives.push_back(payload_buffer);
            next = peel(buffer.data(), payload_buffer.size(), next);
            peeled_count = positives.size() + negatives.size();
        }
    }
    peeled_count = positives.size() + negatives.size();
    if (peeled_count >= peeling_max_iterations) return INFINITE_LOOP;
    auto count_hist = [&](std::array<std::size_t, 4>& hist) {
        for (uint8_t c : counts) {
            for (uint8_t i = 0; i < 4; ++i) {
                ++hist[c & 0x03];
                c >>= 2;
            }
        }
    };
    std::array<std::size_t, 4> histogram = {0,0,0,0};
    count_hist(histogram);
    if (histogram.at(1) or histogram.at(3)) return UNPEELABLE;
    if (histogram.at(2)) return ASYMMETRIC; // full symmetric difference but wrong positive/negative sets depending on the difference order
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

void IBLT::dump_contents() noexcept
{
    std::size_t chunk_size = number_of_buckets / reps;
    bool negative;
    for (std::size_t i = 0; i < number_of_buckets; ++i) {
        fprintf(stderr, "%u | ", get_count_at(i));
        dump_byte_vec(stderr, &hp_buckets.at(i*bucket_size), bucket_size);
        fprintf(stderr, "--> ");
        if (is_peelable(i, negative)) fprintf(stderr, "peelable");
        else fprintf(stderr, "unpeelable");
        fprintf(stderr, "\n");
        if ((i+1) % chunk_size == 0) fprintf(stderr, "\n");
    }
}

//------------------------------------------------------------------------------------------------

std::pair<std::size_t, decltype(IBLT::hasher)::hash_type> IBLT::hash_kmer(uint8_t const* kmer, uint32_t row) const noexcept
{
    auto [hval, cval] = hasher(kmer, payload_buffer.size(), row);
    std::size_t chunk_size = number_of_buckets / reps;
    std::size_t index = (hval % chunk_size) + row * chunk_size;
    return std::make_pair(index, cval);
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
    // std::cerr << "[update count at] ";
    std::size_t byte_idx = idx / 4;
    std::size_t bit_couple_idx = idx % 4;
    uint8_t shift = (3 - bit_couple_idx) * 2;
    uint8_t c = (counts.at(byte_idx) >> shift) & 0x03;
    // std::cerr << "shift = " << uint32_t(shift) << ", [" << byte_idx << "][" << bit_couple_idx << "] = " << uint32_t(counts.at(byte_idx));
    counts[byte_idx] &= ~(0x03 << shift); // clear the 2 bits of the counter
    // std::cerr << " (c: " << uint32_t(c);
    c = op(c);
    c &= 0x03;
    // std::cerr << " --> " << uint32_t(c) << ")";
    counts[byte_idx] |= c << shift;
    // std::cerr << " --> " << uint32_t(counts.at(byte_idx)) << std::endl;
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

void IBLT::xor_at(std::size_t idx, uint8_t const * const kmer, uint64_t hrc) // hrc = hash redundancy code
{
    const std::size_t bucket_byte_start = idx * bucket_size;
    std::size_t payload_shift = hrc_bit_size / 8;

    assert(kmer);
    assert(payload_shift <= sizeof(decltype(hrc)));

    // std::cerr << "[XOR at] ";
    // std::cerr << "idx = " << idx << ", bucket start = " << bucket_byte_start << ", payload shift = " << payload_shift << "\n";
    // std::cerr << "kmer = \n";
    // print_byte_vec(std::cerr, kmer, kmp::ceil_size(2*klen, 8));
    // std::cerr << "\nbucket before addition: \n";
    // print_byte_vec(std::cerr, &hp_buckets[bucket_byte_start], bucket_size);

    uint8_t *const hrc_ptr = reinterpret_cast<uint8_t*>(&hrc);
    // XOR hash value of hash check to bucket (left justified) except for the last byte (if any)
    for (std::size_t i = 0; i < payload_shift; ++i) hp_buckets[i + bucket_byte_start] ^= hrc_ptr[i]; 
    // std::cerr << "\nhash check xored\n";
    if (prefix_len) { // if the hash check is not perfectly aligned, mask and XOR the last byte.
        hrc_ptr[payload_shift] &= mask;
        hp_buckets[payload_shift + bucket_byte_start] ^= hrc_ptr[payload_shift];
        if (bucket_size - ceil_size(2*klen,8) == hash_redundancy_code_buffer.size()) { // hash check and payload share no bytes
            ++payload_shift;
        }
    }
    // XOR k-mer to bucket (right-justified)
    for (std::size_t i = 0; i < payload_buffer.size(); ++i) {
        auto bucket_byte_index = bucket_byte_start + payload_shift + i;
        if (bucket_byte_index >= hp_buckets.size()) {
            throw std::runtime_error(std::string("[IBLT] Overflow error : ") + std::to_string(bucket_byte_index) + " = " + 
            std::to_string(bucket_byte_start) + " + " + 
            std::to_string(payload_shift) + " + " +
            std::to_string(i) + 
            " (prefix length = " + std::to_string(prefix_len) + "," + 
            " hrc byte size = " + std::to_string(hash_redundancy_code_buffer.size()) + ")");
        }
        hp_buckets[bucket_byte_index] ^= kmer[i];
    }
    
    // std::cerr << "bucket after addition: \n";
    // print_byte_vec(std::cerr, &hp_buckets[bucket_byte_start], bucket_size);
    // std::cerr << std::endl;
}

std::size_t IBLT::unpack_at(std::size_t idx) noexcept
{
    assert(number_of_buckets % reps == 0);
    assert(idx < number_of_buckets);
    const std::size_t bucket_byte_start = idx * bucket_size;
    std::size_t payload_shift = hrc_bit_size / 8;

    // std::cerr << "[Unpack at]\n";
    // std::cerr << "\tbucket at " << idx << ":\n";
    // print_byte_vec(std::cerr, &hp_buckets[bucket_byte_start], bucket_size);
    // std::cerr << "\n";

    auto start_itr = hp_buckets.cbegin() + bucket_byte_start;
    hash_redundancy_code_buffer.assign(start_itr, start_itr + hash_redundancy_code_buffer.size());

    bool aligned = prefix_len and (bucket_size - ceil_size(2*klen,8) == hash_redundancy_code_buffer.size());
    if (not aligned) hash_redundancy_code_buffer.back() &= mask; // next: start copying payload from same end byte of hrc
    else ++payload_shift; // otherwise, if everything was perfectly aligned, start from next byte

    // std::cerr << "\thash check buffer:\n";
    // print_byte_vec(std::cerr, hash_redundancy_code_buffer.data(), hash_redundancy_code_buffer.size());
    // std::cerr << "\n";

    assert(bucket_size > payload_shift);
    assert(bucket_size - payload_shift == payload_buffer.size());
    payload_buffer.assign(start_itr + payload_shift, start_itr + bucket_size);
    if (not aligned) payload_buffer.front() &= ~mask;

    // std::cerr << "\tpayload buffer: \n";
    // print_byte_vec(std::cerr, payload_buffer.data(), payload_buffer.size());
    // std::cerr << "\n";

    return idx / (number_of_buckets / reps);
}

bool IBLT::is_peelable(std::size_t index, bool& negative) noexcept
{
    assert(index < number_of_buckets);
    negative = false;
    auto c = get_count_at(index);
    if (c == 0 or c == 2) return false;
    std::size_t row = unpack_at(index); // fill buffers with bucket values
    auto [idx, cval] = hash_kmer(payload_buffer.data(), row);
    if (index != idx) return false;
    auto cval_ptr = reinterpret_cast<uint8_t*>(&cval);
    cval_ptr[hash_redundancy_code_buffer.size()-1] &= mask;
    if (std::equal(hash_redundancy_code_buffer.cbegin(), hash_redundancy_code_buffer.cend(), cval_ptr)) {
        if (c == 3) negative = true;
        return true;
    }
    return false;
}

std::vector<uint8_t> IBLT::get_payload(std::size_t idx) noexcept
{
    unpack_at(idx);
    return payload_buffer;
}

/*
 * This is basically a modified version of remove. New peelable buckets are added to the stack as soon as they are found.
 * Note how peel use the right operation depending on the sign of the bucket to be peeled (operation op) whereas remove always decrements a counter.
 */
IBLT::peel_bucket_state_t IBLT::peel(uint8_t const * const kmer, std::size_t kmer_byte_size, peel_bucket_state_t origin_bucket)
{
    assert(kmer_byte_size == payload_buffer.size());
    // std::cerr << "origin = " << origin_bucket.get_index() << " [";
    uint8_t hit = 0;
    bool pushed = false;
    bool negative;
    peel_bucket_state_t next;
    for (uint8_t i = 0; i < reps; ++i) {
        auto [idx, cval] = hash_kmer(kmer, i);
        if (idx == origin_bucket.get_index()) ++hit;
        // std::cerr << idx << " (" << idx % (number_of_buckets / reps) << ") ";
        update_count_at(idx, origin_bucket.get_update_operation());
        xor_at(idx, kmer, cval);
        if (not pushed and is_peelable(idx, negative)) {
            next = peel_bucket_state_t(idx, negative);
            pushed = true; // This is necessary because all peelable elements coming from the same k-mer are then updated later
        }
    }
    // std::cerr << "]\n";
    assert(hit == 1);
    return next;
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
    os << 
    "k = " << uint32_t(obj.klen) << "\n" <<
    "r = " << uint32_t(obj.reps) << "\n" <<
    "epsilon = " << obj.eps << "\n" <<
    "seed = " << obj.pseed << "\n" <<
    "hash bit-width = " << obj.hrc_bit_size << "\n" <<
    "overlapping byte prefix length = " << obj.prefix_len << "\n" <<
    "mask for overlapping byte = " << uint32_t(obj.mask) << "\n" <<
    "buckets = " << obj.number_of_buckets << "\n" <<
    "bucket size = " << obj.bucket_size;
    return os;
}

std::ostream& operator<<(std::ostream& os, IBLT::failure_t const& err)
{
    switch(err) {
        case IBLT::UNPEELABLE:
            os << "Unpeelable sketch";
            break;
        case IBLT::INFINITE_LOOP:
            os << "Infinite peeling";
            break;
        case IBLT::ASYMMETRIC:
            os << "Asymmetric";
            break;
        default:
            throw std::runtime_error("[Formatting] Unrecognized IBLT error code");
    }
    return os;
}

bool operator==(IBLT const& iblt1, IBLT const& iblt2)
{
    // (1) hash is hard-coded
    // (2) as well as ck_table
    bool params_check = 
        iblt1.klen == iblt2.klen and
        iblt1.reps == iblt2.reps and
        iblt1.eps == iblt2.eps and
        iblt1.max_diff == iblt2.max_diff and
        iblt1.pseed == iblt2.pseed and
        iblt1.hrc_bit_size == iblt2.hrc_bit_size and
        iblt1.prefix_len == iblt2.prefix_len and
        iblt1.mask == iblt2.mask and
        iblt1.number_of_inserted_items == iblt2.number_of_inserted_items and
        iblt1.number_of_buckets == iblt2.number_of_buckets and
        iblt1.bucket_size == iblt2.bucket_size and
        iblt1.counts.size() == iblt2.counts.size() and
        iblt1.hp_buckets.size() == iblt2.hp_buckets.size();
    bool counts_check = true;
    if (params_check)
        for (std::size_t i = 0; i < iblt1.counts.size() and counts_check; ++i) 
            if (iblt1.counts.at(i) != iblt2.counts.at(i)) counts_check = false;
    bool buckets_check = true;
    if (params_check)
        for (std::size_t i = 0; i < iblt1.hp_buckets.size() and buckets_check; ++i) 
            if (iblt1.hp_buckets.at(i) != iblt2.hp_buckets.at(i)) buckets_check = false;
    return params_check and counts_check and buckets_check;
}

bool operator!=(IBLT const& iblt1, IBLT const& iblt2)
{
    return not (iblt1 == iblt2);
}

}
