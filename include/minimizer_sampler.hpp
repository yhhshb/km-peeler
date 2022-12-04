#ifndef MINIMIZER_SAMPLER_HPP
#define MINIMIZER_SAMPLER_HPP

#include <vector>
#include <string>

namespace sampler {

template <class Iterator, typename HashFunction>
class minimizer_sampler
{
    public:
        class const_iterator : public std::iterator<std::forward_iterator_tag, typename Iterator::value_type>
        {
            public:
                typedef typename Iterator::value_type value_type;

                const_iterator(minimizer_sampler const& sampler, Iterator const& start);
                const value_type& operator*() const;
                const value_type& operator++();
                value_type operator++(int);

            private:
                struct mm_pair {
                    value_type mer;
                    typename std::result_of<HashFunction>::type hash_value;
                };

                minimizer_sampler const& parent_sampler;
                Iterator itr_start;
                std::vector<mm_pair> window;
                std::size_t widx, minpos;
                void init();
                void reset_window();
                void find_new_min();
                friend bool operator==(const_iterator const& a, const_iterator const& b);
                friend bool operator!=(const_iterator const& a, const_iterator const& b);
        };
        minimizer_sampler(Iterator const& start, Iterator const& stop, uint16_t window_size);
        const_iterator cbegin() const;
        const_iterator cend() const;
        uint16_t get_w() const;

    private:
        Iterator const itr_start;
        Iterator const itr_stop;
        uint16_t w;
};

template <class Iterator, typename HashFunction>
minimizer_sampler<Iterator, HashFunction>::minimizer_sampler(Iterator const& start, Iterator const& stop, uint16_t window_size)
    : itr_start(start), itr_stop(stop), w(window_size)
{}

template <class Iterator, typename HashFunction>
typename minimizer_sampler<Iterator, HashFunction>::const_iterator minimizer_sampler<Iterator, HashFunction>::cbegin() const
{
    return const_iterator(*this, itr_start);
}

template <class Iterator, typename HashFunction>
typename minimizer_sampler<Iterator, HashFunction>::const_iterator minimizer_sampler<Iterator, HashFunction>::cend() const
{
    return const_iterator(*this, itr_stop);
}

template <class Iterator, typename HashFunction>
uint16_t minimizer_sampler<Iterator, HashFunction>::get_w() const
{
    return w;
}

template <class Iterator, typename HashFunction>
minimizer_sampler<Iterator, HashFunction>::const_iterator::const_iterator(minimizer_sampler const& sampler, Iterator const& start) 
    : parent_sampler(sampler), itr_start(start)
{
    init();
}

template <class Iterator, typename HashFunction>
typename minimizer_sampler<Iterator, HashFunction>::const_iterator::value_type const& minimizer_sampler<Iterator, HashFunction>::const_iterator::operator*() const
{
    return window[minpos];
}

template <class Iterator, typename HashFunction>
typename minimizer_sampler<Iterator, HashFunction>::const_iterator::value_type const& minimizer_sampler<Iterator, HashFunction>::const_iterator::operator++()
{
    auto item = *itr_start++;
    if (not item) {
        reset_window();
    } else {
        bool outside = (widx == minpos);
        window[widx] = {item, HashFunction(item)};
        if (outside) find_new_min();
        else if (window[minpos].hash_value > window[widx].hash_value) minpos = widx;
        widx = (widx + 1) % parent_sampler.w;
    }
    return window[minpos];
}

template <class Iterator, typename HashFunction>
typename minimizer_sampler<Iterator, HashFunction>::const_iterator::value_type minimizer_sampler<Iterator, HashFunction>::const_iterator::operator++(int)
{
    auto current = window[minpos];
    operator++();
    return current;
}

template <class Iterator, typename HashFunction>
void minimizer_sampler<Iterator, HashFunction>::const_iterator::init()
{
    window.resize(parent_sampler.w);
    reset_window();
}

template <class Iterator, typename HashFunction>
void minimizer_sampler<Iterator, HashFunction>::const_iterator::reset_window()
{
    minpos = 0;
    for (std::size_t widx = 0; itr_start != itr_stop and widx < parent_sampler.w; ++widx) {
        auto item = *itr_start++;
        if (not item) {
            widx = 0;
            minpos = 0;
        } else {
            window[widx] = {*item, HashFunction(*item)};
            if (window[widx].hash_value > window[minpos]) minpos = widx;
        }
    }
}

template <class Iterator, typename HashFunction>
void minimizer_sampler<Iterator, HashFunction>::const_iterator::find_new_min()
{
    for (std::size_t i = widx, minpos = widx; i < window.size(); ++i) {
        if (window[minpos].hash_value > window[i].hash_value) minpos = i;
    }
    for (std::size_t i = 0; i < widx; ++i) {
        if (window[minpos].hash_value > window[i].hash_value) minpos = i;
    }
}

template <class Iterator, typename HashFunction>
bool operator==(typename minimizer_sampler<Iterator, HashFunction>::const_iterator const& a, 
                typename minimizer_sampler<Iterator, HashFunction>::const_iterator const& b)
{
    bool same_range = (a.parent_view.itr_start == b.parent_view.itr_start and a.parent_view.itr_end == b.parent_view.itr_end);
    bool same_start = a.itr_start == b.itr_start;
    return same_range and same_start;
}

template <class Iterator, typename HashFunction>
bool operator!=(typename minimizer_sampler<Iterator, HashFunction>::const_iterator const& a, 
                typename minimizer_sampler<Iterator, HashFunction>::const_iterator const& b)
{
    return not (a == b);
}

} // namespace sampler

#endif // MINIMIZER_SAMPLER_HPP