#ifndef CORRECT_HPP
#define CORRECT_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_correct();
int correct_main(const argparse::ArgumentParser& parser);

} // namespace kmp

#endif // CORRECT_HPP