#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include "../include/correct.hpp"
#include "../include/compile_constants.hpp"
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/ordered_unique_sampler.hpp"
#include "../include/IBLT.hpp"
#include "../include/utils.hpp"
#include "../bundled/biolib/include/kmer_view.hpp"
#include <unordered_set>

KSEQ_INIT(gzFile, gzread)

namespace kmp {

int correct_main(const argparse::ArgumentParser& args)
{
    std::size_t loaded_bytes;
    std::string m_filename = args.get<std::string>("extended-minuend");
    std::string mc_filename = args.get<std::string>("correction-minuend");
    std::string s_filename = args.get<std::string>("extended-subtrahend");
    std::string sc_filename = args.get<std::string>("correction-subtrahend");
    uint64_t max_ram = args.get<uint64_t>("--max-ram");
    max_ram *= 1000000000ULL;
    std::string tmp_dir = args.get<std::string>("--tmp-dir");
    bool verbose = args.get<bool>("--verbose");

    IBLT mned = load(m_filename, loaded_bytes);
    if (verbose) std::cerr << "[Info] minuend: loaded " << loaded_bytes << " bytes\n";
    {
        IBLT sbhd = load(s_filename, loaded_bytes);
        if (verbose) std::cerr << "[Info] subtrahend: loaded " << loaded_bytes << " bytes\n";
        mned.subtract(sbhd);
    }
    IBLT cm = load(mc_filename, loaded_bytes);
    if (verbose) std::cerr << "[Info] correction minuend: loaded " << loaded_bytes << " bytes\n";
    {
        IBLT cs = load(sc_filename, loaded_bytes);
        if (verbose) std::cerr << "[Info] correction subtrahend: loaded " << loaded_bytes << " bytes\n";
        cm.subtract(cs);
    }

    std::size_t extended_k = mned.get_k();
    std::size_t k = cm.get_k();
    if (verbose) std::cerr << "[debug] extented k = " << extended_k << ", k = " << k << "\n";
    assert(extended_k >= k);
    kmer_t mask = (2 * k != sizeof(mask) * 8) ? ((kmer_t(1ULL) << (2 * k)) - 1) : std::numeric_limits<decltype(mask)>::max();
    kmp::little2big(reinterpret_cast<uint8_t*>(&mask), sizeof(mask));
    if (verbose) {
        print_byte_vec(std::cerr, reinterpret_cast<uint8_t*>(&mask), sizeof(mask));
        std::cerr << "\n";
    }

    std::vector<std::vector<uint8_t>> p, n;
    IBLT::failure_t err = mned.list(p, n);
    if (verbose) std::cerr << "[debug] extended syncmers: p.size() = " << p.size() << ", n.size() = " << n.size() << std::endl;

    auto get_kmer_set = [&extended_k, &k, &mask](std::vector<std::vector<uint8_t>>& extended_syncmers, emem::external_memory_vector<kmer_t>& kset) {
        for (auto& essb : extended_syncmers) {
            std::size_t shift = sizeof(kmer_t) - essb.size();
            kmer_t ext_syncmer = 0;
            uint8_t* ext_syncmer_ptr = reinterpret_cast<uint8_t*>(&ext_syncmer);
            memcpy(ext_syncmer_ptr + shift, essb.data(), essb.size());
            for (std::size_t i = 0; i < extended_k - k; ++i) {
                kmer_t kmer = ext_syncmer & mask;
                kset.push_back(kmer);
                // cm.remove(reinterpret_cast<uint8_t*>(&kmer), sizeof(kmer)); // remove and insert take care of padding bits at the beginning
                kmp::little2big(ext_syncmer_ptr, sizeof(ext_syncmer));
                ext_syncmer >>= 2;
                kmp::little2big(ext_syncmer_ptr, sizeof(ext_syncmer));
            }
        }
    };
    emem::external_memory_vector<kmer_t> symmetric_difference(max_ram, tmp_dir, "esyncdiff");
    if (not err) {
        emem::external_memory_vector<kmer_t> positive_kmer_set(max_ram, tmp_dir, "positive_kmers");
        get_kmer_set(p, positive_kmer_set);
        sampler::ordered_unique_sampler pusample(positive_kmer_set.cbegin(), positive_kmer_set.cend());
        for (auto kmer : pusample) cm.remove(reinterpret_cast<uint8_t*>(&kmer), sizeof(kmer));

        emem::external_memory_vector<kmer_t> negative_kmer_set(max_ram, tmp_dir, "negative_kmers");
        get_kmer_set(n, negative_kmer_set);
        sampler::ordered_unique_sampler nusample(negative_kmer_set.cbegin(), negative_kmer_set.cend());
        for (auto kmer : nusample) cm.insert(reinterpret_cast<uint8_t*>(&kmer), sizeof(kmer));

        std::set_symmetric_difference(pusample.cbegin(), pusample.cend(), nusample.cbegin(), nusample.cend(), std::back_inserter(symmetric_difference));
    } else {
        std::cerr << "Unable to retrieve extended syncmers\n";
        return err;
    }
    
    // assert(std::adjacent_find(symmetric_difference.cbegin(), symmetric_difference.cend()) == symmetric_difference.cend());

    p.clear();
    n.clear();
    err = cm.list(p, n);
    if (verbose) std::cerr << "[debug] false positive k-mers: p.size() = " << p.size() << ", n.size() = " << n.size() << "\n";
    emem::external_memory_vector<kmer_t> false_positives(max_ram, tmp_dir, "false_positives");
    if (not err) {
        kmer_t km;
        uint8_t* ext_syncmer_ptr = reinterpret_cast<uint8_t*>(&km);
        for (auto& kmer : p) {
            memcpy(ext_syncmer_ptr + sizeof(kmer_t) - kmer.size(), kmer.data(), kmer.size());
            // kmer_t km = *reinterpret_cast<kmer_t*>(kmer.data());
            false_positives.push_back(km & mask);
        }
        for (auto& kmer : n) {
            memcpy(ext_syncmer_ptr + sizeof(kmer_t) - kmer.size(), kmer.data(), kmer.size());
            // kmer_t km = *reinterpret_cast<kmer_t*>(kmer.data());
            false_positives.push_back(km & mask);
        }
    } else {
        std::cerr << "Unable to retrieve false positives: " << err << "\n";
        return err;
    }
    emem::external_memory_vector<kmer_t> true_sym_diff(max_ram, tmp_dir, "true_diff");
    sampler::ordered_unique_sampler unique_diff(symmetric_difference.cbegin(), symmetric_difference.cend());

    // for (auto itr = unique_diff.cbegin(); itr != unique_diff.cend(); ++itr) std::cerr << *itr << "\n";
    // std::cerr << "[debug] unique_diff.size() = " << *unique_diff.size() << "\n";

    // for (auto itr = false_positives.cbegin(); itr != false_positives.cend(); ++itr) std::cerr << *itr << "\n";

    std::set_difference(unique_diff.cbegin(), unique_diff.cend(),
                        false_positives.cbegin(), false_positives.cend(),
                        std::back_inserter(true_sym_diff));
    if (args.is_used("--check")) {
        if (verbose) {
            std::cerr << "[debug] approximated k-mer symmetric difference size = " << *unique_diff.size() << "\n";
            std::cerr << "[debug] false positives (k-mers) size = " << false_positives.size() << "\n";
            std::cerr << "[debug] true k-mer symmetric difference size = " << true_sym_diff.size() << "\n";
        }
        {
            emem::external_memory_vector<kmer_t> subset_check(max_ram, tmp_dir, "subset_check");
            std::set_union(unique_diff.cbegin(), unique_diff.cend(),
                        false_positives.cbegin(), false_positives.cend(),
                        std::back_inserter(subset_check));
            if (subset_check.size() != *unique_diff.size()) {
                throw std::runtime_error("[Check] The false positives are not a subset of the extended syncmer difference");
            }
        }
        {
            emem::external_memory_vector<kmer_t, false> buffer(max_ram, tmp_dir, "buffer");
            {
                gzFile fp;
                kseq_t* seq;
                auto files = args.get<std::vector<std::string>>("--check");
                std::array<emem::external_memory_vector<kmer_t>, 2> kmer_vectors = {
                    emem::external_memory_vector<kmer_t> (max_ram, tmp_dir, "kmer_set1"),
                    emem::external_memory_vector<kmer_t> (max_ram, tmp_dir, "kmer_set2")
                };
                for (std::size_t i = 0; i < files.size(); ++i) {
                    if ((fp = gzopen(files[i].c_str(), "r")) == NULL) throw std::runtime_error("[Check] Unable to open the input file " + files[i] + "\n");
                    seq = kseq_init(fp);
                    while (kseq_read(seq) >= 0) {
                        wrapper::kmer_view<kmer_t> view(seq->seq.s, seq->seq.l, k);
                        for (auto itr = view.cbegin(); itr != view.cend(); ++itr) {
                            if (*itr) kmer_vectors[i].push_back(**itr);
                        }
                    }
                    if (seq) kseq_destroy(seq);
                    gzclose(fp);
                }
                sampler::ordered_unique_sampler kmer_set1(kmer_vectors[0].cbegin(), kmer_vectors[0].cend());
                sampler::ordered_unique_sampler kmer_set2(kmer_vectors[1].cbegin(), kmer_vectors[1].cend());
                std::set_symmetric_difference(kmer_set1.cbegin(), kmer_set2.cend(), 
                                      kmer_set2.cbegin(), kmer_set2.cend(), 
                                      std::back_inserter(buffer));
            }
            emem::external_memory_vector<kmer_t> exact_sym_diff(max_ram, tmp_dir, "exact_diff");
            for (auto itr = buffer.cbegin(); itr != buffer.cend(); ++itr) {
                auto kmer = *itr;
                uint8_t* vptr = reinterpret_cast<uint8_t*>(&kmer);
                kmp::little2big(vptr, sizeof(kmer));
                exact_sym_diff.push_back(*reinterpret_cast<kmer_t*>(vptr));
            }
            if (exact_sym_diff.size() != true_sym_diff.size()) {
                if (verbose) {
                    std::cerr << "[debug] exact symmetric difference size (ground truth) = " << exact_sym_diff.size() << "\n";
                    std::cerr << "[debug] computed symmetric difference size = " << true_sym_diff.size() << "\n";
                }
                throw std::runtime_error("[Check] sizes do not match");
            }
            auto eitr = exact_sym_diff.cbegin();
            auto titr = true_sym_diff.cbegin();
            for (std::size_t i = 0; i < exact_sym_diff.size(); ++i) {
                if (*eitr != *titr) throw std::runtime_error("[Check] spurious elements found");
                ++eitr;
                ++titr;
            }
        }
        std::cerr << "[Check] Everything OK\n";
    }
    if (args.is_used("--list")) {
        std::string list_filename = args.get<std::string>("--list");
        auto cout_buffer_save = std::cout.rdbuf();
        std::ofstream myostrm;
        if (list_filename != ".") { // redirect cout
            myostrm.open(list_filename);
            std::cout.rdbuf(myostrm.rdbuf());
        }

        for (auto itr = true_sym_diff.cbegin(); itr != true_sym_diff.cend(); ++itr) {
            std::cout << *itr << "\n";
        }

        std::cout.rdbuf(cout_buffer_save); // restore cout
    }
    return 0;
}

argparse::ArgumentParser get_parser_correct()
{
    argparse::ArgumentParser parser("correct");
    parser.add_description("Check correctness of the combination between extended syncmers and whole IBLTs");
    parser.add_argument("extended-minuend")
        .help("extended IBLT from which to subtract");
    parser.add_argument("correction-minuend")
        .help("correction IBLT of the minuend (containing all k-mers)");
    parser.add_argument("extended-subtrahend")
        .help("extended IBLT to subtract");
    parser.add_argument("correction-subtrahend")
        .help("correction IBLT of the subtrahend (containing all k-mers)");
    parser.add_argument("-l", "--list")
        .help("List IBLT to file (use '.' for stdout)");
    parser.add_argument("-m", "--max-ram")
        .help("RAM limit")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(6));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-c", "--check")
        .help("Original fasta files used to build the two IBLTs (only syncmers are supported)")
        .nargs(2);
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

} // namespace kmp