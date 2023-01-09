#ifndef ORDERED_UNIQUE_SAMPLER_HPP
#define ORDERED_UNIQUE_SAMPLER_HPP

#include <vector>

namespace sampler {

template <class Iterator>
class ordered_unique_sampler
{
    public:
        class const_iterator
        {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type   = std::ptrdiff_t;
                using value_type        = typename Iterator::value_type;
                using pointer           = value_type*;
                using reference         = value_type&;

                const_iterator(ordered_unique_sampler const& sampler, Iterator const& start);
                value_type const& operator*() const;
                const_iterator const& operator++();
                const_iterator operator++(int);

            private:
                // struct mm_pair {
                //     value_type mer;
                //     typename std::result_of<HashFunction>::type hash_value;
                // };

                ordered_unique_sampler const& parent_sampler;
                Iterator itr_start;
                value_type buffer;
                friend bool operator==(const_iterator const& a, const_iterator const& b);
                friend bool operator!=(const_iterator const& a, const_iterator const& b);
        };

        ordered_unique_sampler(Iterator const& start, Iterator const& stop);
        const_iterator cbegin() const;
        const_iterator cend() const;

    private:
        Iterator const itr_start;
        Iterator const itr_stop;
};

template <class Iterator>
ordered_unique_sampler<Iterator>::ordered_unique_sampler(Iterator const& start, Iterator const& stop)
    : itr_start(start), itr_stop(stop)
{}

template <class Iterator>
typename ordered_unique_sampler<Iterator>::const_iterator ordered_unique_sampler<Iterator>::cbegin() const
{
    return const_iterator(*this, itr_start);
}

template <class Iterator>
typename ordered_unique_sampler<Iterator>::const_iterator ordered_unique_sampler<Iterator>::cend() const
{
    return const_iterator(*this, itr_stop);
}

template <class Iterator>
ordered_unique_sampler<Iterator>::const_iterator::const_iterator(ordered_unique_sampler const& sampler, Iterator const& start) 
    : parent_sampler(sampler), itr_start(start)
{}

template <class Iterator>
typename ordered_unique_sampler<Iterator>::const_iterator::value_type const& ordered_unique_sampler<Iterator>::const_iterator::operator*() const
{
    return *itr_start;
}

template <class Iterator>
typename ordered_unique_sampler<Iterator>::const_iterator const& ordered_unique_sampler<Iterator>::const_iterator::operator++()
{
    auto prev = *itr_start;
    while(prev == *itr_start && itr_start != parent_sampler.itr_stop) ++itr_start;
    return *this;
}

template <class Iterator>
typename ordered_unique_sampler<Iterator>::const_iterator ordered_unique_sampler<Iterator>::const_iterator::operator++(int)
{
    auto current = *this;
    operator++();
    return current;
}

} // namespace sampler

#endif // ORDERED_UNIQUE_SAMPLER_HPP