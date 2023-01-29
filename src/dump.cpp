#include "../include/dump.hpp"
#include "../include/IBLT.hpp"

namespace kmp {

int dump_main(const argparse::ArgumentParser& args)
{
    std::size_t loaded_bytes;
    IBLT sketch = load(args.get<std::string>("sketch"), loaded_bytes);
    std::cerr << "[Info] minuend: loaded " << loaded_bytes << " bytes\n";
    sketch.dump_contents();
    std::cout << "\n";
    return 0;
}

argparse::ArgumentParser get_parser_dump()
{
    argparse::ArgumentParser parser("dump");
    parser.add_description("Dump IBLT in a human-readable way");
    parser.add_argument("sketch")
        .help("the name of the IBLT to print")
        .required();
    return parser;
}

}