// #include <cstdio>
#include "../include/diff.hpp"
#include "../include/io.hpp"
#include "../include/IBLT.hpp"
#include "../include/utils.hpp"

namespace kmp {

int diff_main(const argparse::ArgumentParser& args)
{
    std::size_t loaded_bytes, m_size, s_size;
    std::string m_filename = args.get<std::string>("minuend");
    std::string s_filename = args.get<std::string>("subtrahend");
    IBLT mned = load(m_filename, loaded_bytes);
    m_size = mned.size();
    std::cerr << "[Info] minuend: loaded " << loaded_bytes << " bytes\n";
    // mned.print_config(std::cerr);
    { // memory management of sbhd
        IBLT sbhd = load(s_filename, loaded_bytes);
        s_size = sbhd.size();
        std::cerr << "[Info] subtrahend: loaded " << loaded_bytes << " bytes\n";
        // sbhd.print_config(std::cerr);
        mned.subtract(sbhd);
    }
    
    if (args.is_used("--output")) io::store(mned, args.get<std::string>("--output"));
    if (args.is_used("--jaccard")) {
        std::size_t m_diff_size, s_diff_size;
        IBLT::failure_t err;
        auto cout_buffer_save = std::cout.rdbuf();
        std::ofstream myostrm;
        if (args.get<std::string>("--jaccard") != ".") { // redirect cout
            myostrm.open(args.get<std::string>("--jaccard"));
            std::cout.rdbuf(myostrm.rdbuf());
        }
        std::cout << m_filename << "," << s_filename << ",";
        if (not (err = mned.list(m_diff_size, s_diff_size))) {
            auto intersection_size = m_size - m_diff_size;
            assert(intersection_size == s_size - s_diff_size);
            double j = double(intersection_size) / double(m_size + s_size - intersection_size);
            std::cout << j;
        } else {
            std::cout << err;
        }
        std::cout << "\n";
        std::cout.rdbuf(cout_buffer_save); // restore cout
    }
    if (args.is_used("--list")) {
        std::vector<std::vector<uint8_t>> p, n;
        IBLT::failure_t err;
        auto cout_buffer_save = std::cout.rdbuf();
        std::ofstream myostrm;
        if (args.get<std::string>("--list") != ".") { // redirect cout
            myostrm.open(args.get<std::string>("--list"));
            std::cout.rdbuf(myostrm.rdbuf());
        }
        if (not (err = mned.list(p, n))) {
            for (auto& v : p) std::cout << vec2kmer(v, mned.get_k()) << "\n";
            for (auto& v : n) std::cout << vec2kmer(v, mned.get_k()) << "\n";
        } else {
            std::cout << err;
        }
        std::cout.rdbuf(cout_buffer_save); // restore cout
    }
    return 0;
}

argparse::ArgumentParser get_parser_diff()
{
    argparse::ArgumentParser parser("diff");
    parser.add_argument("minuend")
        .help("IBLT from which to subtract");
    parser.add_argument("subtrahend")
        .help("IBLT to subtract");
    parser.add_argument("-o", "--output")
        .help("Save difference IBLT to disk");
    parser.add_argument("-j", "--jaccard")
        .help("Compute Jaccard similarity from IBLT difference. Save result as CSV line (minuend, subtrahend, jaccard). Use '.' for stdout")
        .default_value("");
    parser.add_argument("-l", "--list")
        .help("List IBLT to file. Use '.' for stdout");
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

}