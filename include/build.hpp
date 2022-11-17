#ifndef MAINS_HPP
#define MAINS_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_build();
int build_main(const argparse::ArgumentParser& parser);

} // namespace kmp

#endif // MAINS_HPP