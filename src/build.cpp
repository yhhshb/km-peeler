#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include <exception>
#include <iostream>

#include "../include/compile_constants.hpp"
#include "../include/build.hpp"
#include "../include/external_memory_vector.hpp"
#include "../include/kmer_view.hpp"
#include "../include/syncmer_sampler.hpp"
#include "../include/ordered_unique_sampler.hpp"
#include "../include/hash.hpp"
#include "../include/io.hpp"

#include "../include/IBLT.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kmp {

int build_main(const argparse::ArgumentParser& args)
{
    gzFile fp;
    kseq_t* seq;

    std::string input_filename = args.get<std::string>("input_file");
    std::string output_filename = args.get<std::string>("output_filename");
    uint64_t n = args.get<uint64_t>("-n");
    uint8_t k = args.get<uint8_t>("k");
    uint8_t z = args.get<uint8_t>("z");
    uint8_t offset = args.get<uint8_t>("--offset");
    uint8_t r = args.get<uint8_t>("--repetitions");
    float epsilon = args.get<float>("--epsilon");
    std::string tmp_dir = args.get<std::string>("--tmp-dir");
    uint64_t seed = args.get<uint64_t>("--seed");
    uint64_t max_ram = args.get<uint64_t>("--max-ram");
    uint64_t max_ram_bytes = max_ram * 1000000000ULL;
    bool verbose = args.get<bool>("--verbose");

    std::size_t k_max = (sizeof(kmer_t) * 8 / 2);
    if (k > k_max) throw std::invalid_argument("k cannot be greater than " + std::to_string(k_max));
    if (z > k) throw std::invalid_argument("z cannot be greater than k");
    if (r > 7 or r < 3) throw std::invalid_argument("r must be an integer value in [3, 7]");
    if (epsilon > 1 or epsilon < 0) throw std::invalid_argument("epsilon must be a floating point number in [0, 1]");

    if (verbose) std::cerr << "Part 1: file reading and info gathering\n";

    emem::external_memory_vector<kmer_t> kmer_vector(max_ram_bytes, tmp_dir, "kmers");
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) throw std::runtime_error("Unable to open the input file " + input_filename + "\n");
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        wrapper::kmer_view<kmer_t> view(seq->seq.s, seq->seq.l);
        sampler::syncmer_sampler<wrapper::kmer_view<kmer_t>::const_iterator, hash::mm_pos_extractor> syncmers(view.cbegin(), view.cend(), z, offset);
        for (auto itr = syncmers.cbegin(); itr != syncmers.cend(); ++itr) {
            // if (*itr) kmer_vector.push_back(**itr); // kmer_view iterators return std::optionals
            kmer_vector.push_back(*itr);
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);

    IBLT iblt(k, r, epsilon, n, seed);
    sampler::ordered_unique_sampler unique_kmers(kmer_vector.cbegin(), kmer_vector.cend());
    for (auto itr = unique_kmers.cbegin(); itr != unique_kmers.cend(); ++itr) {
        iblt.insert(reinterpret_cast<uint8_t const * const>(&(*itr)), sizeof(*itr));
    }

    std::cerr << "Written IBLT of " << io::store(iblt, output_filename) << " Bytes\n";

    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_argument("input_file")
        .help("fastx file in input (gzipped or not)");
    parser.add_argument("output_file")
        .help("output file name");
    parser.add_argument("-n")
        .help("Expected number of differences in the IBLT")
        .scan<'i', uint64_t>();
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'d', uint8_t>()
        .default_value(0);
    parser.add_argument("-z")
        .help("syncmer parameter")
        .scan<'d', uint8_t>()
        .default_value(0);
    parser.add_argument("-f", "--offset")
        .help("z-mer offset inside syncmers")
        .scan<'d', uint8_t>()
        .default_value(0);
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
