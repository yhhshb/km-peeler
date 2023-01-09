#include <iostream>
#include <argparse/argparse.hpp>
#include "../include/build.hpp"

using namespace kmp;

void print_help()
{
    std::cerr << "Available sub-commands are:\n";
    std::cerr << "\tbuild\n";
    std::cerr << "\tTODO\n";
}

int main(int argc, char *argv[])
{
    auto build_parser = get_parser_build();
    // auto diff_parser = get_parser_diff();
    // auto list_parser = get_list_build();

    argparse::ArgumentParser program("kmp");
    program.add_subparser(build_parser);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(program);
    // else if (strcmp(argv[1], "diff") == 0) return list_main(argc - 1, &argv[1]);
    // else if (strcmp(argv[1], "list") == 0) return diff_main(argc - 1, &argv[1]);
    else throw std::runtime_error("[kmp.cpp] This should never happen");
    return 0;
}

// cmake .. -D CMAKE_BUILD_TYPE=Debug -D KMP_USE_SANITIZERS=On
// make -p