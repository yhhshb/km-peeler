#include <tuple>
#include "../include/bias.hpp"
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/ordered_unique_sampler.hpp"
#include "../bundled/biolib/include/kmer_view.hpp"
#include "../bundled/biolib/include/minimizer_sampler.hpp"
#include "../bundled/biolib/include/syncmer_sampler.hpp"
#include "../bundled/biolib/include/hash_sampler.hpp"
#include "../bundled/biolib/include/sequence_generator.hpp"
#include "../bundled/biolib/include/sequence_mutator.hpp"
#include "../bundled/biolib/include/jaccard.hpp"
#include "../bundled/biolib/include/hash.hpp"

namespace kmp {

typedef sequence::generator::nucleic<sequence::generator::packed::uniform> ug_t;

std::tuple<double, double, double, double> test_sampling(ug_t&, uint64_t, uint64_t, uint8_t, uint8_t, double, bool, std::string, uint64_t, bool);

int bias_main(const argparse::ArgumentParser& parser)
{
    auto l = parser.get<uint64_t>("--length");
    auto k = parser.get<uint16_t>("-k");
    auto z = parser.get<uint16_t>("-z");
    auto max_rate = parser.get<double>("--max-mutation-rate");
    auto nsteps = parser.get<uint64_t>("--steps");
    auto trials = parser.get<uint64_t>("--trials");
    auto canonical = parser.get<bool>("--canonical");
    auto tmp_dir = parser.get<std::string>("--tmp-dir");
    auto seed = parser.get<uint64_t>("--seed");
    auto max_ram = parser.get<uint64_t>("--max-ram");
    auto verbose = parser.get<bool>("--verbose");

    if (k > 32) throw std::invalid_argument("0 <= k <= 32");
    if (z > k) throw std::invalid_argument("z <= k");
    if (max_rate < 0 or max_rate > 1) throw std::invalid_argument("0 <= --max-rate <= 1");
    
    sequence::generator::packed::uniform uniform_sequence_generator(seed);
    sequence::generator::nucleic generator(uniform_sequence_generator);

    std::random_device dev;
    std::mt19937 rng(dev());
    // std::uniform_int_distribution<std::mt19937::result_type> dist6(1,6);

    double increment = max_rate / nsteps;
    for (std::size_t i = 0; i < nsteps + 1; ++i) {
        double mut_rate = increment * i;
        double exc, hbc, mbc, sbc, err_hb, err_mm, err_sb, sqerr_hb, sqerr_mm, sqerr_sb;
        exc = hbc = mbc = sbc = err_hb = err_mm = err_sb = sqerr_hb = sqerr_mm = sqerr_sb = 0;
        for (std::size_t j = 0; j < trials; ++j) {
            auto [exact, hash_based, mm_based, sync_based] = test_sampling(generator, rng(), l, k, z, mut_rate, canonical, tmp_dir, max_ram, verbose);
            exc += exact;
            hbc += hash_based;
            mbc += mm_based;
            sbc += sync_based;
            err_hb = hash_based - exact;
            err_mm = mm_based - exact;
            err_sb = sync_based - exact;
            sqerr_hb += err_hb * err_hb;
            sqerr_mm += err_mm * err_mm;
            sqerr_sb += err_sb * err_sb;
        }
        std::cout << 
            exc / trials << "," << 
            hbc / trials << "," << 
            mbc / trials << "," << 
            sbc / trials << "," << 
            sqerr_hb / trials << "," << 
            sqerr_mm / trials << "," << 
            sqerr_sb /trials << "\n";
    }
    return 0;
}

std::tuple<double, double, double, double> test_sampling(ug_t& generator, uint64_t str_gen_seed, uint64_t l, uint8_t k, uint8_t z, double mut_rate, bool canonical, std::string tmp_dir, uint64_t max_ram, bool verbose)
{
    double exact, hash, mm, sync, sampling_rate_A, sampling_rate_B;
    std::string original = generator.get_sequence(l);
    wrapper::sequence_mutator mutator_view(original.begin(), original.end(), mut_rate, 0, 0);
    std::string mutated;
    decltype(mutator_view)::const_iterator::report_t stats;
    for (auto itr = mutator_view.cbegin(str_gen_seed); itr != mutator_view.cend(); ++itr) {
        if (*itr) mutated.push_back(**itr);
        stats = itr.get_report();
    }
    uint64_t emem_size = max_ram / 2 == 0 ? 1000 : max_ram / 2 * 1000;
    std::size_t kmer_size_A, kmer_size_B;
    { // memory collection
        emem::external_memory_vector<uint64_t> A(emem_size, tmp_dir, "SetAfull");
        emem::external_memory_vector<uint64_t> B(emem_size, tmp_dir, "SetBfull");
        wrapper::kmer_view<uint64_t> A_view(original, k, canonical);
        wrapper::kmer_view<uint64_t> B_view(mutated, k, canonical);
        for (auto itr = A_view.cbegin(); itr != A_view.cend(); ++itr) if (*itr) A.push_back(**itr);
        for (auto itr = B_view.cbegin(); itr != B_view.cend(); ++itr) if (*itr) B.push_back(**itr);

        sampler::ordered_unique_sampler unique_kmers_A(A.cbegin(), A.cend());
        sampler::ordered_unique_sampler unique_kmers_B(B.cbegin(), B.cend());
        auto [num, den, A_size, B_size] = algorithm::jaccard(unique_kmers_A.cbegin(), unique_kmers_A.cend(), unique_kmers_B.cbegin(), unique_kmers_B.cend());
        exact = static_cast<double>(num) / den;

        { // syncmer jaccard
            uint16_t offset = (k - z + 1) / 2;
            hash::minimizer_position_extractor mmp_extractor(k, z);
            sampler::syncmer_sampler A_sync(unique_kmers_A.cbegin(), unique_kmers_A.cend(), mmp_extractor, offset, offset);
            sampler::syncmer_sampler B_sync(unique_kmers_B.cbegin(), unique_kmers_B.cend(), mmp_extractor, offset, offset);
            auto [n, d, A_ssize, B_ssize] = algorithm::jaccard(A_sync.cbegin(), A_sync.cend(), B_sync.cbegin(), B_sync.cend());
            sync = static_cast<double>(n) / d;
            sampling_rate_A = static_cast<double>(A_ssize) / A_size;
            sampling_rate_B = static_cast<double>(B_ssize) / B_size;
            if (verbose) {
                std::cerr << "Original syncmer sampling rate: " << sampling_rate_A << "\n";
                std::cerr << "Mutated  syncmer sampling rate: " << sampling_rate_B << "\n"; 
            }
        }
        { // hash based
            sampler::hash_sampler A_hash(unique_kmers_A.cbegin(), unique_kmers_A.cend(), hash::hash64(), 0, sampling_rate_A);
            sampler::hash_sampler B_hash(unique_kmers_B.cbegin(), unique_kmers_B.cend(), hash::hash64(), 0, sampling_rate_A); // use the same sampling rate
            auto [n, d, A_hsize, B_hsize] = algorithm::jaccard(A_hash.cbegin(), A_hash.cend(), B_hash.cbegin(), B_hash.cend());
            hash = static_cast<double>(n) / d;
            if (verbose) {
                std::cerr << "Original hash sampling rate: " << static_cast<double>(A_hsize) / A_size << "\n";
                std::cerr << "Mutated  hash sampling rate: " << static_cast<double>(B_hsize) / B_size << "\n"; 
            }
        }
        kmer_size_A = A_size; // carry-over to minimizer section
        kmer_size_B = B_size;
    } // end of memory collection
    { // minimizer based
        emem::external_memory_vector<uint64_t> A(emem_size, tmp_dir, "SetAmm");
        emem::external_memory_vector<uint64_t> B(emem_size, tmp_dir, "SetBmm");
        uint16_t window_size = static_cast<uint16_t>(2.0 / sampling_rate_A - 1);
        wrapper::kmer_view<uint64_t> A_view(original, k, canonical);
        wrapper::kmer_view<uint64_t> B_view(mutated, k, canonical);
        sampler::minimizer_sampler A_mm(A_view.cbegin(), A_view.cend(), hash::hash64(), 0, window_size);
        sampler::minimizer_sampler B_mm(B_view.cbegin(), B_view.cend(), hash::hash64(), 0, window_size);
        for (auto itr = A_mm.cbegin(); itr != A_mm.cend(); ++itr) if (*itr) A.push_back(**itr);
        for (auto itr = B_mm.cbegin(); itr != B_mm.cend(); ++itr) if (*itr) B.push_back(**itr);
        sampler::ordered_unique_sampler unique_mm_A(A.cbegin(), A.cend());
        sampler::ordered_unique_sampler unique_mm_B(B.cbegin(), B.cend());
        auto [n, d, A_msize, B_msize] = algorithm::jaccard(unique_mm_A.cbegin(), unique_mm_A.cend(), unique_mm_B.cbegin(), unique_mm_B.cend());
        mm = static_cast<double>(n) / d;
        if (verbose) {
            std::cerr << "Original hash sampling rate: " << static_cast<double>(A_msize) / kmer_size_A << "\n";
            std::cerr << "Mutated  hash sampling rate: " << static_cast<double>(B_msize) / kmer_size_B << "\n"; 
        }
    }
    return std::make_tuple(exact, hash, mm, sync);
}

argparse::ArgumentParser get_parser_bias()
{
    argparse::ArgumentParser parser("bias");
    parser.add_description("Test unbiased behaviour of (open) syncmers (minimum z-mer in the middle)");
    // parser.add_argument("output_filename")
    //     .help("output file name");
    parser.add_argument("-l", "--length")
        .help("sequence length")
        .scan<'i', uint64_t>()
        .required();
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'u', uint16_t>()
        .required();
    parser.add_argument("-z")
        .help("syncmer parameter")
        .scan<'u', uint16_t>()
        .required();
    parser.add_argument("-x", "--max-mutation-rate")
        .help("maximum mutation rate to test")
        .scan<'g', double>()
        .required();
    parser.add_argument("-n", "--steps")
        .help("number of mutation rates to test (form 0 to --max-mutation-rate)")
        .scan<'i', uint64_t>()
        .required();
    parser.add_argument("-t", "--trials")
        .help("number of trials for each point")
        .scan<'i', uint64_t>()
        .required();
    parser.add_argument("-c", "--canonical")
        .help("Canonical k-mers")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-s", "--seed")
        .help("random seed")
        .scan<'i', uint64_t>()
        .default_value(uint64_t(42));
    parser.add_argument("-m", "--max-ram")
        .help("RAM limit in MB")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(1000));
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

} // namespace kmp