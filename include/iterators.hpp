#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <vector>
#include <string>

namespace kmp {

template <std::size_t KmerLen>
class kmer_view
{
    public:
        struct kmer_t {
            uint8_t kmer[(KmerLen + 8 - 1 / 8)];
        };

        class const_iterator : public std::iterator<std::forward_iterator_tag, kmer_t>
        {
            public:
                typedef kmer_t value_type;

                const_iterator();
                const_iterator(kmer_view const* view);
                const kmer_t& operator*();
                const kmer_t& operator++();
                kmer_t operator++(int);
            
            private:
                kmer_view const* parent_view;
                std::size_t char_idx;
                uint8_t mask;
                kmer_t buffer;
                std::size_t bases_since_last_break;
                void find_first_good_kmer();
                friend bool operator==(const_iterator const& a, const_iterator const& b);
                friend bool operator!=(const_iterator const& a, const_iterator const& b);
        };
    
        kmer_view(char const* contig, std::size_t contig_len);
        kmer_view(std::string const& contig);
        const_iterator cbegin() const;
        const_iterator cend() const;
        std::size_t get_k() const;
        std::size_t get_kmer_bit_width() const;

    private:
        char const* seq;
        std::size_t slen;
        std::size_t klen;
        std::size_t bit_width;
};

template <std::size_t KmerLen>
kmer_view<KmerLen>::kmer_view(char const* contig, std::size_t contig_len) 
: seq(contig), slen(contig_len), klen(KmerLen), bit_width((KmerLen + 8 - 1 / 8))
{
    // klen = static_cast<decltype(klen)>(bit_width + 8 - 1) / 8;
}

template <class Iterator>
class minimizer_sampler
{
    public:
        class const_iterator : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
        {
            public:
                const_iterator();
                const_iterator(minimizer_sampler const* sampler);
                const kmer_t& operator*();
                const kmer_t& operator++();
                kmer_t operator++(int);
        };
        minimizer_sampler(Iterator const& km_itr, std::size_t window_size);
        const_iterator begin() const;
        const_iterator end() const;
        uint16_t get_w() const;

    private:
        Iterator const& itr;
        std::size_t w;
};

template <class Iterator>
class syncmer_sampler
{
    public:
        class const_iterator : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
        {
            public:
                typedef typename Iterator::value_type value_type;
                const_iterator();
                const_iterator(syncmer_sampler const* sampler);
                const value_type& operator*();
                const value_type& operator++();
                value_type operator++(int);
        };

        syncmer_sampler(Iterator const& km_itr, uint16_t substring_size, uint16_t offset);
        const_iterator begin() const;
        const_iterator end() const;
        uint16_t get_w() const;

    private:
        Iterator const& itr;
        uint16_t z;
        uint16_t offset;
};

// template <class Iterator>
// syncmer_sampler<Iterator>::syncmer_sampler(Iterator const& km_itr, uint16_t substring_size, uint16_t offset) 
//     : itr(km_itr), z(substring_size), offset(offset)
// {

// }

// template <class Iterator>
// syncmer_sampler<Iterator>::const_iterator syncmer_sampler<Iterator>::begin() const
// {
    
// }

// template <class Iterator>
// syncmer_sampler<Iterator>::const_iterator syncmer_sampler<Iterator>::end() const
// {

// }

// template <class Iterator>
// uint16_t syncmer_sampler<Iterator>::get_w() const
// {

// }

// template <class Iterator>
// syncmer_sampler<Iterator>::const_iterator::const_iterator()
// {

// }

// template <class Iterator>
// syncmer_sampler<Iterator>::const_iterator::const_iterator(syncmer_sampler const* sampler)
// {

// }

// template <class Iterator>
// const Iterator::value_type& syncmer_sampler<Iterator>::const_iterator::operator*()
// {

// }

// template <class Iterator>
// const Iterator::value_type& syncmer_sampler<Iterator>::const_iterator::operator++()
// {

// }

// template <class Iterator>
// typename Iterator::value_type syncmer_sampler<Iterator>::const_iterator::operator++(int)
// {

// }

} // namespace kmp

#endif // SAMPLER_HPP