#ifndef IBLT_HPP
#define IBLT_HPP

#include <array>
#include <vector>
#include <string>
#include <stack>
#include <functional>
#include "../include/hash.hpp"

namespace kmp {

class IBLT
{
    public:
        enum failure_t {NONE, UNPEELABLE, INFINITE_LOOP, ASYMMETRIC};
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
        void dump_contents() noexcept;
        template <class Visitor> void visit(Visitor& visitor);

    private:
        class peel_bucket_state_t {
            public:
                peel_bucket_state_t() : index(std::numeric_limits<decltype(index)>::max()), negative(true) {}
                peel_bucket_state_t(std::size_t idx, bool neg) : index(idx), negative(neg) {}
                std::size_t get_index() const noexcept {return index;}
                bool is_valid() const noexcept {return index != std::numeric_limits<decltype(index)>::max();}
                bool is_negative() const noexcept {return negative;}
                std::function<uint8_t(uint8_t)> get_update_operation() const noexcept {return negative ? [](uint8_t c){return c + 1;} : [](uint8_t c){return c - 1;};}
            private:
                std::size_t index;
                bool negative;
        };
        const std::array<float, 8> ck_table = {0, 0, 0, 1.222, 1.295, 1.425, 1.570, 1.721};
        //const removed since we need loading
        uint8_t klen;
        uint8_t reps;
        double eps;
        uint64_t max_diff;
        uint64_t pseed;
        uint16_t hrc_bit_size;
        std::size_t prefix_len;
        uint8_t mask;
        // end const block
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
        std::pair<std::size_t, decltype(hasher)::hash_type> hash_kmer(uint8_t const* kmer, uint32_t row) const noexcept;
        uint8_t get_count_at(std::size_t idx) const noexcept;
        void update_count_at(std::size_t idx, std::function<uint8_t(uint8_t)>) noexcept;
        void inc_count_at(std::size_t idx) noexcept;
        void dec_count_at(std::size_t idx) noexcept;
        void xor_at(std::size_t idx, uint8_t const * const kmer, uint64_t header);
        std::size_t unpack_at(std::size_t idx) noexcept;
        bool is_peelable(std::size_t idx, bool& negative) noexcept;
        std::vector<uint8_t> get_payload(std::size_t idx) noexcept;
        peel_bucket_state_t peel(uint8_t const * const kmer, std::size_t kmer_byte_size, peel_bucket_state_t origin_bucket);

        friend IBLT load(std::string filename, std::size_t& byte_size);
        friend std::ostream& operator<<(std::ostream& os, IBLT const& obj);
        friend bool operator==(IBLT const& iblt1, IBLT const& iblt2);
};

template <class Visitor>
void IBLT::visit(Visitor& visitor)
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
bool operator==(IBLT const& iblt1, IBLT const& iblt2);
bool operator!=(IBLT const& iblt1, IBLT const& iblt2);

} // namespace kmp

#endif // IBLT_HPP