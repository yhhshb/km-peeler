#ifndef BIAS_HPP
#define BIAS_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_bias();
int bias_main(const argparse::ArgumentParser& parser);

} // namespace kmp

#endif // BIAS_HPP