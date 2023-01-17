#ifndef BUILD_HPP
#define BUILD_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_build();
int build_main(const argparse::ArgumentParser& parser);

} // namespace kmp

#endif // BUILD_HPP