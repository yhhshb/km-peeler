#ifndef IBLT_HPP
#define IBLT_HPP

#include <array>
#include <vector>
#include <string>

namespace kmp {

class IBLT
{
    public:

        IBLT(uint8_t k, uint8_t r, double epsilon, uint64_t n, uint64_t seed);
        virtual void insert(std::vector<uint8_t> const& kmer) noexcept;
        virtual void remove(std::vector<uint8_t> const& kmer) noexcept;
        void subtract(const IBLT& v);
        std::vector<std::vector<uint8_t>> list();
        ~IBLT();

    private:
        const std::array<float, 8> ck_table = {0, 0, 0, 1.222, 1.295, 1.425, 1.570, 1.721};
        uint8_t klen;
        uint8_t reps;
        double eps;
        uint64_t max_diff;
        uint64_t pseed;
        uint64_t number_of_buckets;
        uint16_t hf_size;
        uint32_t bucket_size;
        std::vector<uint8_t> counts;
        std::vector<uint8_t> hp_buckets;

        void update_count_at(std::size_t idx, std::function<uint8_t(uint8_t)>) noexcept;
        void inc_count_at(std::size_t idx) noexcept;
        void dec_count_at(std::size_t idx) noexcept;
        void xor_at(std::size_t idx, std::vector<uint8_t> const& kmer) noexcept;
        std::size_t find_peelable_bucket() const;
        std::vector<uint8_t> get_payload(std::size_t idx) const noexcept;

        friend std::ostream& operator<<(std::ostream& os, const IBLT& obj);
};

std::ostream& operator<<(std::ostream& os, const IBLT& obj);

}

#endif // IBLT_HPP