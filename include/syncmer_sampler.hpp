#ifndef SYNCMER_SAMPLER_HPP
#define SYNCMER_SAMPLER_HPP

#include <vector>

namespace sampler {

template <class Iterator, typename PropertyExtractor>
class syncmer_sampler
{
    public:
        class const_iterator : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
        {
            public:
                typedef typename Iterator::value_type value_type;
                
                const_iterator(syncmer_sampler const& sampler, Iterator const& start);
                const value_type& operator*() const;
                const value_type& operator++();
                value_type operator++(int);

            private:
                syncmer_sampler const& parent_sampler;
                Iterator itr_start;
        };

        syncmer_sampler(Iterator const& start, Iterator const& stop, uint16_t substring_size, uint16_t offset);
        const_iterator begin() const;
        const_iterator end() const;
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
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator syncmer_sampler<Iterator, PropertyExtractor>::begin() const
{
    return const_iterator(*this, itr_start);
}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator syncmer_sampler<Iterator, PropertyExtractor>::end() const
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
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::value_type const& syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::operator++()
{
    while (PropertyExtractor(*itr_start, z) != soffset) ++itr_start;
}

template <class Iterator, typename PropertyExtractor>
typename syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::value_type syncmer_sampler<Iterator, PropertyExtractor>::const_iterator::operator++(int)
{
    auto current = *itr_start;
    operator++();
    return current;
}

} // namespace sampler

#endif // SYNCMER_SAMPLER_HPP