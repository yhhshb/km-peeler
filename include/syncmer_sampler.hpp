#ifndef SYNCMER_SAMPLER_HPP
#define SYNCMER_SAMPLER_HPP

#include <vector>

namespace sampler {

template <class Iterator, typename PropertyExtractor>
class syncmer_sampler
{
    public:
        class const_iterator
        {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type   = std::ptrdiff_t;
                using value_type        = typename Iterator::value_type::value_type;
                using pointer           = value_type*;
                using reference         = value_type&;
                
                const_iterator(syncmer_sampler const& sampler, Iterator const& start);
                value_type const& operator*() const;
                const_iterator const& operator++();
                const_iterator operator++(int);

            private:
                syncmer_sampler const& parent_sampler;
                Iterator itr_start;
                friend bool operator==(const_iterator const& a, const_iterator const& b);
                friend bool operator!=(const_iterator const& a, const_iterator const& b);
        };

        syncmer_sampler(Iterator const& start, Iterator const& stop, uint16_t substring_size, uint16_t offset);
        const_iterator cbegin() const;
        const_iterator cend() const;
        uint16_t get_substring_size() const;

    private:
        Iterator const itr_start;
        Iterator const itr_stop;
        uint16_t z;
        uint16_t soffset;
};

template <class Iterator, typename PropertyExtractor>
syncmer_sampler<Iterator, PropertyExtractor>::syncmer_sampler(Iterator const& start, Iterator const& stop, uint16_t substring_size, uint16_t offset) 
    : itr_start(start), itr_stop(stop), z(substring_size), soffset(offset)
{}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator syncmer_sampler<Iterator, PropertyExtractor>::cbegin() const
{
    return const_iterator(*this, itr_start);
}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator syncmer_sampler<Iterator, PropertyExtractor>::cend() const
{
    return const_iterator(*this, itr_stop);
}

template <class Iterator, typename PropertyExtractor>
uint16_t syncmer_sampler<Iterator, PropertyExtractor>::get_substring_size() const
{
    return z;
}

template <class Iterator, typename PropertyExtractor>
syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::const_iterator(syncmer_sampler const& sampler, Iterator const& start)
    : parent_sampler(sampler), itr_start(start)
{}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::value_type const& syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::operator*() const
{
    return *itr_start;
}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator const& syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::operator++()
{
    while (PropertyExtractor(*itr_start, parent_sampler.z) != parent_sampler.soffset) ++itr_start;
    return *this;
}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::operator++(int)
{
    auto current = *itr_start;
    operator++();
    return current;
}

template <class Iterator, typename HashFunction>
bool operator==(typename syncmer_sampler<Iterator, HashFunction>::const_iterator const& a, 
                typename syncmer_sampler<Iterator, HashFunction>::const_iterator const& b)
{
    bool same_range = (a.parent_view.itr_start == b.parent_view.itr_start and a.parent_view.itr_end == b.parent_view.itr_end);
    bool same_start = a.itr_start == b.itr_start;
    return same_range and same_start;
}

template <class Iterator, typename HashFunction>
bool operator!=(typename syncmer_sampler<Iterator, HashFunction>::const_iterator const& a, 
                typename syncmer_sampler<Iterator, HashFunction>::const_iterator const& b)
{
    return not (a == b);
}

} // namespace sampler

#endif // SYNCMER_SAMPLER_HPP