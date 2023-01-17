#ifndef DIFF_HPP
#define DIFF_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_diff();
int diff_main(const argparse::ArgumentParser& parser);

} // namespace kmp

#endif // DIFF_HPP