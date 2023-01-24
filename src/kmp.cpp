#include <iostream>
#include <argparse/argparse.hpp>
#include "../include/build.hpp"
#include "../include/diff.hpp"
#include "../include/correct.hpp"

using namespace kmp;

int main(int argc, char *argv[])
{
    auto build_parser = get_parser_build();
    auto diff_parser = get_parser_diff();
    auto correct_parser = get_parser_correct();

    argparse::ArgumentParser program("kmp");
    program.add_subparser(build_parser);
    program.add_subparser(diff_parser);
    program.add_subparser(correct_parser);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else if (program.is_subcommand_used(diff_parser)) return diff_main(diff_parser);
    else if (program.is_subcommand_used(correct_parser)) return correct_main(program);
    else std::cerr << program << std::endl;
    return 0;
}

// mkdir debug_build
// cd debug_build
// cmake .. -D CMAKE_BUILD_TYPE=Debug -D KMP_USE_SANITIZERS=On
// make