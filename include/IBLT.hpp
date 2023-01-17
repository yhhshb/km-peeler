#ifndef IBLT_HPP
#define IBLT_HPP

#include <array>
#include <vector>
#include <string>
#include <stack>
#include "../include/hash.hpp"

namespace kmp {

class IBLT
{
    public:
        enum failure_t {NONE, UNPEELABLE, INFINITE_LOOP};
        IBLT(uint8_t k, uint8_t r, double epsilon, uint64_t n, uint64_t seed);
        IBLT (IBLT &) = default;
        IBLT (IBLT &&) = default;
        void insert(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept;
        void remove(uint8_t const * const kmer, std::size_t kmer_byte_size) noexcept;
        void subtract(IBLT const& other);
        failure_t list(std::vector<std::vector<uint8_t>>& positives, std::vector<std::vector<uint8_t>>& negatives);
        failure_t list(std::size_t& positive_size, std::size_t& negative_size);
        std::size_t size() const noexcept;
        uint8_t get_k() const noexcept;
        template <class Visitor> void visit(Visitor const& visitor);

    private:
        typedef std::stack<std::size_t, std::vector<std::size_t>> idx_stack_t;
        const std::array<float, 8> ck_table = {0, 0, 0, 1.222, 1.295, 1.425, 1.570, 1.721};
        const uint8_t klen;
        const uint8_t reps;
        const double eps;
        const uint64_t max_diff;
        const uint64_t pseed;
        const uint16_t hrc_bit_size;
        const std::size_t prefix_len;
        const uint8_t mask;
        std::size_t number_of_inserted_items;
        uint64_t number_of_buckets;
        uint32_t bucket_size;
        std::vector<uint8_t> counts;
        std::vector<uint8_t> hp_buckets;
        std::vector<uint8_t> hash_redundancy_code_buffer; // not saved to disk
        std::vector<uint8_t> payload_buffer; // not saved to disk
        hash::double_hash64 hasher;

        IBLT();
        void resize_buffers();
        uint8_t get_count_at(std::size_t idx) const noexcept;
        void update_count_at(std::size_t idx, std::function<uint8_t(uint8_t)>) noexcept;
        void inc_count_at(std::size_t idx) noexcept;
        void dec_count_at(std::size_t idx) noexcept;
        void xor_at(std::size_t idx, uint8_t const * const kmer, uint64_t header) noexcept;
        std::size_t unpack_at(std::size_t idx) noexcept;
        bool is_peelable(std::size_t idx) noexcept;
        std::size_t find_peelable_bucket() noexcept;
        std::vector<uint8_t> get_payload(std::size_t idx) noexcept;
        void peel(uint8_t const * const kmer, std::size_t kmer_byte_size, std::size_t origin_bucket, idx_stack_t& idxs);

        friend IBLT load(std::string filename, std::size_t& byte_size);
        friend std::ostream& operator<<(std::ostream& os, IBLT const& obj);
};

template <class Visitor>
void IBLT::visit(Visitor const& visitor)
{
    visitor.apply(klen);
    visitor.apply(reps);
    visitor.apply(eps);
    visitor.apply(max_diff);
    visitor.apply(pseed);
    visitor.apply(hrc_bit_size);
    visitor.apply(prefix_len);
    visitor.apply(mask);
    visitor.apply(number_of_inserted_items);
    visitor.apply(number_of_buckets);
    visitor.apply(bucket_size);
    visitor.apply(counts);
    visitor.apply(hp_buckets);
}

IBLT load(std::string filename, std::size_t& byte_size);
std::ostream& operator<<(std::ostream& os, IBLT const& obj);
std::ostream& operator<<(std::ostream& os, IBLT::failure_t const& obj);

} // namespace kmp

#endif // IBLT_HPP