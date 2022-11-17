#include "iterators.hpp"
#include "utils.hpp"
#include <array>

namespace kmp {

const std::array<uint8_t, 256> seq_nt4_table = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

kmer_view::kmer_view(char const* contig, std::size_t contig_len, uint16_t k) : seq(contig), len(contig_len), klen(k)
{
    //
}

kmer_view::const_iterator kmer_view::begin() const
{
    return const_iterator(this);
}

kmer_view::const_iterator kmer_view::end() const
{
    return const_iterator();
}

uint16_t kmer_view::get_k() const
{
    return klen;
}

kmer_view::const_iterator::const_iterator() : parent_view(nullptr), char_idx(0), bases_since_last_break(0) {}

kmer_view::const_iterator::const_iterator(kmer_view const* view) : parent_view(view), char_idx(0), bases_since_last_break(0)
{
    uint16_t buf_size = parent_view->klen / 4 + (parent_view->klen % 4 ? 1 : 0);
    uint16_t padding = buf_size * 8 - 2 * parent_view->klen;
    assert(padding <= 6);
    mask = ~((uint8_t(1) << (8 - padding)) - 1);
    buffer.resize(buf_size);
    find_first_good_kmer();
}

const std::vector<uint8_t>& kmer_view::const_iterator::operator*()
{
    return buffer;
}

const std::vector<uint8_t>& kmer_view::const_iterator::operator++()
{
    assert(char_idx <= parent_view->len);
    if (char_idx == parent_view->len) parent_view = nullptr;
    assert(parent_view);
    char c = parent_view->seq[char_idx++];
    if (seq_nt4_table[c] < 4) {
        shift_byte_vec(buffer.data(), buffer.data(), buffer.size(), 2);
        buffer[0] &= mask;
        buffer.back() |= c;
        ++bases_since_last_break;
    } else {
        bases_since_last_break = 0;
        find_first_good_kmer();
    }
    return operator*();
}

std::vector<uint8_t> kmer_view::const_iterator::operator++(int)
{
    auto res = buffer;
    operator++();
    return res;
}

void kmer_view::const_iterator::find_first_good_kmer()
{
    while(char_idx < parent_view->len and bases_since_last_break < parent_view->klen) {
        char c = parent_view->seq[char_idx++];
        if (seq_nt4_table[c] < 4) {
            shift_byte_vec(buffer.data(), buffer.data(), buffer.size(), 2);
            buffer[0] &= mask;
            buffer.back() |= c;
            ++bases_since_last_break;
        } else {
            bases_since_last_break = 0;
        }
    }
}

bool operator==(kmer_view::const_iterator const& a, kmer_view::const_iterator const& b)
{
    if (a.parent_view == nullptr and b.parent_view == nullptr) return true;
    return a.char_idx == b.char_idx;
}

bool operator!=(kmer_view::const_iterator const& a, kmer_view::const_iterator const& b)
{
    return not (a == b);
}

}