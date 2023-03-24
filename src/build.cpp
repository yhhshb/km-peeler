#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include <exception>
#include <iostream>

#include "../include/build.hpp"
#include "../include/compile_constants.hpp"
#include "../include/external_memory_vector.hpp"
#include "../include/kmer_view.hpp"
#include "../include/syncmer_sampler.hpp"
#include "../include/ordered_unique_sampler.hpp"
#include "../include/hash.hpp"
#include "../include/io.hpp"
#include "../include/IBLT.hpp"
#include "../include/utils.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kmp {

int build_main(const argparse::ArgumentParser& args)
{
    gzFile fp;
    kseq_t* seq;

    std::string input_filename = args.get<std::string>("input_filename");
    std::string output_filename = args.get<std::string>("output_filename");
    uint64_t n = args.get<uint64_t>("-n");
    uint8_t k = args.get<uint8_t>("-k");
    uint8_t z = args.get<uint8_t>("-z");
    uint8_t offset1 = args.get<uint8_t>("--offset");
    uint8_t x = args.get<uint8_t>("--extend");
    uint8_t r = args.get<uint8_t>("--repetitions");
    float epsilon = args.get<float>("--epsilon");
    std::string tmp_dir = args.get<std::string>("--tmp-dir");
    uint64_t seed = args.get<uint64_t>("--seed");
    uint64_t max_ram = args.get<uint64_t>("--max-ram");
    uint64_t max_ram_bytes = max_ram * 1000000000ULL;
    bool canonical = args.get<bool>("--canonical");
    bool verbose = args.get<bool>("--verbose");

    std::size_t k_max = (sizeof(kmer_t) * 8 / 2);
    if (k == 0) throw std::invalid_argument("k == 0");
    if ((k + x) > k_max) throw std::invalid_argument("(k + x) > " + std::to_string(k_max));
    if (z > k) throw std::invalid_argument("z > k");
    if (r > 7 or r < 3) throw std::invalid_argument("r must be an integer value in [3, 7]");
    if (epsilon > 1 or epsilon < 0) throw std::invalid_argument("epsilon must be a floating point number in [0, 1]");

    uint8_t offset2 = offset1 ? offset1 : k - z; // if open syncmers offset2 is the same as offset1

    if (verbose) {
        std::cerr << "Part 1: file reading and info gathering\n";
        std::cerr << "\toffsets: (" << uint32_t(offset1) << ", " << uint32_t(offset2) << ")\n";
        std::cerr << "\textention = " << uint32_t(x) << "\n";
    }

    emem::external_memory_vector<kmer_t> kmer_vector(max_ram_bytes, tmp_dir, "kmers");
    hash::minimizer_position_extractor mmp_extractor(k, z);
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) throw std::runtime_error("Unable to open the input file " + input_filename + "\n");
    seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        wrapper::kmer_view<kmer_t> view(seq->seq.s, seq->seq.l, k+x, canonical);
        sampler::syncmer_sampler syncmers(view.cbegin(), view.cend(), mmp_extractor, offset1, offset2);
        for (auto itr = syncmers.cbegin(); itr != syncmers.cend(); ++itr) {
            kmer_vector.push_back(*itr);
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);

    IBLT iblt(k + x, r, epsilon, n, seed);
    if (verbose) {
        std::cerr << "Part 2: found " << kmer_vector.size() << " syncmers\n";
        std::cerr << iblt << "\n\n";
    }
    sampler::ordered_unique_sampler unique_kmers(kmer_vector.cbegin(), kmer_vector.cend());
    std::size_t unique_kmers_size = 0;
    for (auto itr = unique_kmers.cbegin(); itr != unique_kmers.cend(); ++itr) {
        auto val = *itr;
        auto vptr = reinterpret_cast<uint8_t*>(&val);
        kmp::little2big(vptr, sizeof(val));
        iblt.insert(vptr, sizeof(val)); // FIXME hide endian check inside insert
        ++unique_kmers_size;
    }
    if (verbose) std::cerr << "Found " << unique_kmers_size << " unique syncmers\n";
    std::cerr << "Written IBLT of " << io::store(iblt, output_filename) << " Bytes\n";

    if (args.get<bool>("--check")) {
        std::size_t loaded_bytes = 0;
        IBLT loaded = load(output_filename, loaded_bytes);
        std::cerr << "Re-loaded " << loaded_bytes << " bytes\n";
        if (loaded != iblt) std::cerr << "FAILED";
        else std::cerr << "Everything is OK";
        std::cerr << std::endl;
    }
    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_description("Build an Invertible Bloom Lookup Table");
    parser.add_argument("input_filename")
        .help("fastx file in input (gzipped or not)");
    parser.add_argument("output_filename")
        .help("output file name");
    parser.add_argument("-n")
        .help("Expected number of differences in the IBLT")
        .scan<'u', uint64_t>()
        .required();
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'u', uint8_t>()
        .required();
        // .default_value(uint8_t(0)); // NOTE: casting to the correct type here is **VERY** important
    parser.add_argument("-z")
        .help("syncmer parameter")
        .scan<'u', uint8_t>()
        .required();
        // .default_value(uint8_t(0));
    parser.add_argument("-f", "--offset")
        .help("z-mer offset inside syncmers")
        .scan<'u', uint8_t>()
        .default_value(uint8_t(0));
    parser.add_argument("-x", "--extend")
        .help("extension of x bases after each syncmer")
        .scan<'u', uint8_t>()
        .default_value(uint8_t(0));
    parser.add_argument("-c", "--canonical")
        .help("Canonical k-mers")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-r", "--repetitions")
        .help("number of hash functions for the IBLT")
        .scan<'u', uint8_t>()
        .default_value(uint8_t(4));
    parser.add_argument("-e", "--epsilon")
        .help("fraction of space overhead")
        .scan<'g', float>()
        .default_value(float(0));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-s", "--seed")
        .help("random seed")
        .scan<'i', uint64_t>()
        .default_value(uint64_t(42));
    parser.add_argument("-m", "--max-ram")
        .help("RAM limit")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(6));
    parser.add_argument("-C", "--check")
        .help("Check saving/loading consistency")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

} // namespace kmp
