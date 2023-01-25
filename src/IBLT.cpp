#include "../include/IBLT.hpp"

#include <cmath>
#include "../include/io.hpp"
#include "../include/utils.hpp"

#include <iostream>

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
    const std::size_t peeling_max_iterations = 2 * max_diff;
    std::size_t peeled_count = 0;
    long long next;
    // idx_stack_t peelable_idx_stack;
    std::vector<uint8_t> buffer;
    buffer.resize(payload_buffer.size());
    // std::cerr << "[List] ";
    while (
        (next = find_peelable_bucket()) != static_cast<long long>(number_of_buckets) and 
         peeled_count < peeling_max_iterations) // FIXME
    {
        // peelable_idx_stack.push(tmp);
        // std::cerr << "(found peelable bucket at " << tmp << ")\n";
        // while (not peelable_idx_stack.empty()) {
        while (next != static_cast<long long>(number_of_buckets) and peeled_count < peeling_max_iterations) {
            // auto idx = peelable_idx_stack.top();
            // peelable_idx_stack.pop();
            assert(buffer.size() == payload_buffer.size());
            buffer = payload_buffer; // make local copy of kmer since payload_buffer is a shared "global" method variable
            if (next < 0) {
                unpack_at(-next);
                negatives.push_back(payload_buffer);
                next = peel(buffer.data(), payload_buffer.size(), -next, [](uint8_t c){return c + 1;}); //, peelable_idx_stack);
            } else {
                unpack_at(next);
                positives.push_back(payload_buffer);
                next = peel(buffer.data(), payload_buffer.size(), next, [](uint8_t c){return c - 1;}); //, peelable_idx_stack);
            }
            peeled_count = positives.size() + negatives.size();
        }
        // peeled_count = positives.size() + negatives.size();
    }
    peeled_count = positives.size() + negatives.size();
    if (peeled_count >= peeling_max_iterations) return INFINITE_LOOP;
    for (auto c : counts) if (c != 0) return UNPEELABLE;
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

void IBLT::print_config(std::ostream& os) const noexcept
{
    os << 
    "k = " << uint32_t(klen) << "\n" <<
    "r = " << uint32_t(reps) << "\n" <<
    "epsilon = " << eps << "\n" <<
    "seed = " << pseed << "\n" <<
    "hash bit-width = " << hrc_bit_size << "\n" <<
    "overlapping byte prefix length = " << prefix_len << "\n" <<
    "mask for overlapping byte = " << uint32_t(mask) << "\n" <<
    "buckets = " << number_of_buckets << "\n" <<
    "bucket size = " << bucket_size << "\n\n";
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

int IBLT::is_peelable(std::size_t index) noexcept
{
    assert(index < number_of_buckets);
    auto c = get_count_at(index);
    if (c == 0 or c == 2) return 0;
    std::size_t row = unpack_at(index); // fill buffers with bucket values
    auto [idx, cval] = hash_kmer(payload_buffer.data(), row);
    // std::cerr << " [is peelable] row(" << index << ") = " << row << ", index check = " << idx;
    if (index != idx) return 0;
    auto cval_ptr = reinterpret_cast<uint8_t*>(&cval);
    cval_ptr[hash_redundancy_code_buffer.size()-1] &= mask;
    if (std::equal(hash_redundancy_code_buffer.cbegin(), hash_redundancy_code_buffer.cend(), cval_ptr)) 
        return c == 1 ? 1 : -1; // check hrc and recomputed hash value
    return 0;
}

long long IBLT::find_peelable_bucket() noexcept
{
    int sign;
    for (std::size_t i = 0; i < number_of_buckets; ++i) {
        if ((sign = is_peelable(i))) return sign * i;
    }
    return static_cast<long long>(number_of_buckets);
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
long long IBLT::peel(uint8_t const * const kmer, std::size_t kmer_byte_size, std::size_t origin_bucket, std::function<uint8_t(uint8_t)> op) //, idx_stack_t& idxs)
{
    assert(kmer_byte_size == payload_buffer.size());
    // std::cerr << "origin = " << origin_bucket << "[\n";
    uint8_t hit = 0;
    bool pushed = false;
    long long next = static_cast<long long>(number_of_buckets);
    for (uint8_t i = 0; i < reps; ++i) {
        auto [idx, cval] = hash_kmer(kmer, i);
        if (idx == origin_bucket) ++hit;
        // std::cerr << "\tcounts[" << idx << "] = " << uint32_t(get_count_at(idx)) << " -> ";
        update_count_at(idx, op);
        // std::cerr << uint32_t(get_count_at(idx));
        xor_at(idx, kmer, cval);
        int sign;
        if (not pushed and (sign = is_peelable(idx))) {
            // if (sign > 0) std::cerr << "(+) ";
            // else std::cerr << "(-) ";
            // idxs.push(
            next = sign * static_cast<long long>(idx);//);
            pushed = true; // This is necessary because all peelable elements coming from the same k-mer are then updated later
        }
        // std::cerr << "]\n";
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
    // write obj to stream
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
        default:
            throw std::runtime_error("[Formatting] Unrecognized IBLT error code");
    }
    return os;
}

}
