#include <iostream>
#include <argparse/argparse.hpp>
#include "../include/build.hpp"
#include "../include/diff.hpp"
#include "../include/correct.hpp"
#include "../include/bias.hpp"
#include "../include/dump.hpp"

using namespace kmp;

int main(int argc, char *argv[])
{
    auto build_parser = get_parser_build();
    auto diff_parser = get_parser_diff();
    auto correct_parser = get_parser_correct();
    auto bias_parser = get_parser_bias();
    auto dump_parser = get_parser_dump();

    argparse::ArgumentParser program(argv[0]);
    program.add_subparser(build_parser);
    program.add_subparser(diff_parser);
    program.add_subparser(correct_parser);
    program.add_subparser(bias_parser);
    program.add_subparser(dump_parser);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else if (program.is_subcommand_used(diff_parser)) return diff_main(diff_parser);
    else if (program.is_subcommand_used(correct_parser)) return correct_main(correct_parser);
    else if (program.is_subcommand_used(bias_parser)) return bias_main(bias_parser);
    else if (program.is_subcommand_used(dump_parser)) return dump_main(dump_parser);
    else std::cerr << program << std::endl;
    return 0;
}

// mkdir debug_build
// cd debug_build
// cmake .. -D CMAKE_BUILD_TYPE=Debug -D KMP_USE_SANITIZERS=On
// make