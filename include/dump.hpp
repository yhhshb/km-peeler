#ifndef DUMP_HPP
#define DUMP_HPP

#include <argparse/argparse.hpp>

namespace kmp {

argparse::ArgumentParser get_parser_dump();
int dump_main(const argparse::ArgumentParser& parser);

}

#endif