#ifndef KMER_VIEW_HPP
#define KMER_VIEW_HPP

#include <string>
#include <optional>

#include "constants.hpp"

namespace wrapper {

template <typename KmerType>
class kmer_view
{
    public:
        class const_iterator : public std::iterator<std::forward_iterator_tag, KmerType>
        {
            public:
                typedef KmerType value_type;

                const_iterator(kmer_view const* view) noexcept;
                std::optional<KmerType> const& operator*() const noexcept;
                const_iterator const& operator++();
                const_iterator operator++(int);
            
            protected:
                const_iterator(kmer_view const* view, int dummy_end) noexcept;

            private:
                kmer_view const* parent_view;
                std::size_t char_idx;
                KmerType mask;
                std::size_t shift;
                KmerType kmer_buffer[2];
                uint8_t strand;
                std::size_t bases_since_last_break;
                void find_first_good_kmer();
                friend bool operator==(const_iterator const& a, const_iterator const& b);
                friend bool operator!=(const_iterator const& a, const_iterator const& b);
        };
    
        kmer_view(char const* contig, std::size_t contig_len, uint8_t k, bool canonical = false) noexcept;
        kmer_view(std::string const& contig, uint8_t k, bool canonical = false) noexcept;
        const_iterator cbegin() const;
        const_iterator cend() const noexcept;
        uint8_t get_k() const noexcept;

    private:
        char const* seq;
        std::size_t slen;
        uint8_t klen;
        bool canon;
};

//---------------------------------------------------------------------------------------------------------------------------------------------------

template <typename KmerType>
kmer_view<KmerType>::kmer_view(char const* contig, std::size_t contig_len, uint8_t k, bool canonical) noexcept
    : seq(contig), slen(contig_len), klen(k), canon(canonical)
{}

template <typename KmerType>
kmer_view<KmerType>::kmer_view(std::string const& contig, uint8_t k, bool canonical) noexcept
    : seq(contig.c_str()), slen(contig.length()), klen(k), canon(canonical)
{}

template <typename KmerType>
kmer_view<KmerType>::const_iterator kmer_view<KmerType>::cbegin() const
{
    return const_iterator(this);
}

template <typename KmerType>
kmer_view<KmerType>::const_iterator kmer_view<KmerType>::cend() const noexcept
{
    return const_iterator();
}

template <typename KmerType>
uint8_t kmer_view<KmerType>::get_k() const noexcept
{
    return klen;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

template <typename KmerType>
kmer_view<KmerType>::const_iterator::const_iterator(kmer_view const* view, int dummy_end) noexcept
    : parent_view(view), char_idx(view->slen), strand(0), bases_since_last_break(0) 
{}

template <typename KmerType>
kmer_view<KmerType>::const_iterator::const_iterator(kmer_view const* view) noexcept
    : parent_view(view), char_idx(0), strand(0), bases_since_last_break(0)
{
    mask = (KmerType(1) << (2 * parent_view->klen)) - 1;
    find_first_good_kmer();
}

template <typename KmerType>
std::optional<KmerType> const& kmer_view<KmerType>::const_iterator::operator*() const noexcept
{
    if (bases_since_last_break == 0) return std::nullopt;
    return kmer_buffer.at(strand);
}

template <typename KmerType>
kmer_view<KmerType>::const_iterator const& 
kmer_view<KmerType>::const_iterator::operator++()
{
    assert(parent_view);
    assert(char_idx <= parent_view->slen);
    if (bases_since_last_break == 0) {
        find_first_good_kmer();
        return operator*();
    }
    uint8_t c = constants::seq_nt4_table[parent_view->seq[char_idx++]];
    if (char_idx == parent_view->slen) parent_view = nullptr;
    if (c < 4) {
        kmer_buffer[0] = (kmer_buffer[0] << 2 | c) & mask;            /* forward m-mer */
        kmer_buffer[1] = (kmer_buffer[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
        if (parent_view->canon and kmer_buffer[0] != kmer_buffer[1]) strand = kmer_buffer[0] < kmer_buffer[1] ? 0 : 1;  // strand, if symmetric k-mer then use previous strand
        ++bases_since_last_break;
    } else {
        bases_since_last_break = 0;
        // find_first_good_kmer();
    }
    return *this;
}

template <typename KmerType>
kmer_view<KmerType>::const_iterator 
kmer_view<KmerType>::const_iterator::operator++(int)
{
    auto res = *this;
    operator++();
    return res;
}

template <typename KmerType>
void kmer_view<KmerType>::const_iterator::find_first_good_kmer()
{
    while(char_idx < parent_view->slen and bases_since_last_break < parent_view->klen) {
        uint8_t c = constants::seq_nt4_table.at(parent_view->seq[char_idx++]);
        if (c < 4) {
            kmer_buffer[0] = (kmer_buffer[0] << 2 | c) & mask;            /* forward m-mer */
            kmer_buffer[1] = (kmer_buffer[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
            if (parent_view->canon and kmer_buffer[0] != kmer_buffer[1]) strand = kmer_buffer[0] < kmer_buffer[1] ? 0 : 1;
            ++bases_since_last_break;
        } else {
            bases_since_last_break = 0;
        }
    }
}

template <typename KmerType>
bool operator==(typename kmer_view<KmerType>::const_iterator const& a, typename kmer_view<KmerType>::const_iterator const& b)
{
    return a.parent_view == b.parent_view and a.char_idx == b.char_idx;
}

template <typename KmerType>
bool operator!=(typename kmer_view<KmerType>::const_iterator const& a, typename kmer_view<KmerType>::const_iterator const& b)
{
    return not (a == b);
}

} // namespace wrapper

#endif // KMER_VIEW_HPP