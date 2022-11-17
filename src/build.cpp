#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include <exception>
#include <iostream>

#include "../include/build.hpp"
#include "../include/external_memory_vector.hpp"
#include "../include/iterators.hpp"

#include "../include/compile_constants.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kmp {

int build_main(const argparse::ArgumentParser& args)
{
    typedef kmer_view<KLEN>::kmer_t kmer_t;
    gzFile fp;
    kseq_t* seq;

    std::string input_filename = args.get<std::string>("input_file");
    // uint8_t k = args.get<uint8_t>("k");
    std::string tmp_dir = args.get<std::string>("--tmp-dir");
    uint64_t max_ram = args.get<uint64_t>("--max-ram");
    uint64_t max_ram_bytes = max_ram * 1000000000ULL;
    bool verbose = args.get<bool>("--verbose");

    // std::size_t k_max = (sizeof(kmer_t) * 8 / 2);
    // if (k > k_max) throw std::invalid_argument("k cannot be greater than " + std::to_string(k_max));

    if (verbose) std::cerr << "Part 1: file reading and info gathering\n";

    emem::external_memory_vector<kmer_t> kmer_vector(max_ram_bytes, tmp_dir, "kmers");
    std::size_t kmer_width_bits = 0;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) throw std::runtime_error("Unable to open the input file " + input_filename + "\n");
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        kmer_view<KLEN> view(seq->seq.s, seq->seq.l);
        for (auto itr = view.cbegin(); itr != view.cend(); ++itr) {
            kmer_vector.push_back(*itr);
        }
        assert(view.get_k() == KLEN);
        auto prev_bw = kmer_width_bits;
        kmer_width_bits = view.get_kmer_bit_width();
        if (prev_bw and prev_bw != kmer_width_bits) throw std::runtime_error("Not all k-mer view have the same k-mer type");
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);

    for (auto itr = kmer_vector.cbegin(); itr != kmer_vector.cend(); ++itr) {
        // remove duplicates
        // sample
        // insert into IBLT
    }
    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_argument("input_file")
        .help("fastx file in input (gzipped or not)");
    parser.add_argument("output_file")
        .help("output file name");
    // parser.add_argument("-k")
    //     .help("k-mer length")
    //     .scan<'d', uint8_t>()
    //     .default_value(0);
    parser.add_argument("-r", "--repetitions")
        .help("number of hash functions for the IBLT")
        .scan<'d', uint8_t>()
        .default_value(3);
    parser.add_argument("-e", "--epsilon")
        .help("fraction of space overhead")
        .scan<'g', float>()
        .default_value(0);
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(".");
    parser.add_argument("-s", "--seed")
        .help("random seed")
        .scan<'i', uint64_t>()
        .default_value(42);
    parser.add_argument("-m", "--max-memory")
        .help("RAM limit")
        .scan<'d', uint64_t>()
        .default_value(6);
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

} // namespace kmp
