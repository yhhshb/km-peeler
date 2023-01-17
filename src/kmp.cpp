#include <iostream>
#include <argparse/argparse.hpp>
#include "../include/build.hpp"
#include "../include/diff.hpp"

using namespace kmp;

// void print_help()
// {
//     std::cerr << "Available sub-commands are:\n";
//     std::cerr << "\tbuild\n";
//     std::cerr << "\tdiff\n";
// }

int main(int argc, char *argv[])
{
    auto build_parser = get_parser_build();
    auto diff_parser = get_parser_diff();
    // auto list_parser = get_list_build();

    argparse::ArgumentParser program("kmp");
    program.add_subparser(build_parser);
    program.add_subparser(diff_parser);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else if (program.is_subcommand_used(diff_parser)) return diff_main(diff_parser);
    // else if (strcmp(argv[1], "list") == 0) return list_main(program);
    else std::cerr << program << std::endl;
    return 0;
}

// mkdir debug_build
// cd debug_build
// cmake .. -D CMAKE_BUILD_TYPE=Debug -D KMP_USE_SANITIZERS=On
// make