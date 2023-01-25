#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include "../include/correct.hpp"
#include "../include/compile_constants.hpp"
#include "../include/external_memory_vector.hpp"
#include "../include/IBLT.hpp"
#include "../include/utils.hpp"
#include "../include/kmer_view.hpp"

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
    IBLT mned = load(m_filename, loaded_bytes);
    std::cerr << "[Info] minuend: loaded " << loaded_bytes << " bytes\n";
    {
        IBLT sbhd = load(s_filename, loaded_bytes);
        std::cerr << "[Info] subtrahend: loaded " << loaded_bytes << " bytes\n";
        mned.subtract(sbhd);
    }
    IBLT cm = load(mc_filename, loaded_bytes);
    std::cerr << "[Info] correction minuend: loaded " << loaded_bytes << " bytes\n";
    {
        IBLT cs = load(sc_filename, loaded_bytes);
        std::cerr << "[Info] correction subtrahend: loaded " << loaded_bytes << " bytes\n";
        cm.subtract(cs);
    }
    std::vector<std::vector<uint8_t>> p, n;
    IBLT::failure_t err = mned.list(p, n);
    std::size_t extended_k = mned.get_k();
    std::size_t k = cm.get_k();
    assert(extended_k >= k);
    kmer_t mask;
    if (2 * k != sizeof(mask) * 8) mask = (kmer_t(1) << (2 * k)) - 1;
    else mask = std::numeric_limits<decltype(mask)>::max();

    emem::external_memory_vector<kmer_t> symmetric_difference(max_ram, tmp_dir, "esyncdiff");
    if (not err) {
        for (auto& essb : p) {
            kmer_t ext_syncmer = *reinterpret_cast<kmer_t*>(essb.data());
            for (std::size_t i = 0; i < extended_k - k; ++i) {
                kmer_t kmer = ext_syncmer & mask;
                symmetric_difference.push_back(kmer);
                cm.remove(reinterpret_cast<uint8_t*>(&kmer), sizeof(kmer_t));
                ext_syncmer >>= 2;
            }
        }
        for (auto& essb : n) {
            kmer_t ext_syncmer = *reinterpret_cast<kmer_t*>(essb.data());
            for (std::size_t i = 0; i < extended_k - k; ++i) {
                kmer_t kmer = ext_syncmer & mask;
                symmetric_difference.push_back(kmer);
                cm.insert(reinterpret_cast<uint8_t*>(&ext_syncmer), sizeof(kmer_t));
                ext_syncmer >>= 2;
            }
        }
    } else {
        std::cerr << "Unable to retrieve extended syncmers\n";
        return err;
    }

    p.clear();
    n.clear();
    err = cm.list(p, n);
    emem::external_memory_vector<kmer_t> false_positives(max_ram, tmp_dir, "esyncdiff");
    if (not err) {
        for (auto& kmer : p) {
            kmer_t km = *reinterpret_cast<kmer_t*>(kmer.data());
            false_positives.push_back(km);
        }
        for (auto& kmer : n) {
            kmer_t km = *reinterpret_cast<kmer_t*>(kmer.data());
            false_positives.push_back(km);
        }
    } else {
        std::cerr << "Unable to retrieve false positives\n";
        return err;
    }
    emem::external_memory_vector<kmer_t> true_sym_diff(max_ram, tmp_dir, "true_diff");
    std::set_difference(symmetric_difference.cbegin(), symmetric_difference.cend(),
                        false_positives.cbegin(), false_positives.cend(),
                        std::back_inserter(true_sym_diff));
    if (args.is_used("--check")) {
        {
            emem::external_memory_vector<kmer_t> subset_check(max_ram, tmp_dir, "subset_check");
            std::set_union(symmetric_difference.cbegin(), symmetric_difference.cend(),
                        false_positives.cbegin(), false_positives.cend(),
                        std::back_inserter(subset_check));
            if (subset_check.size() != symmetric_difference.size())
                throw std::runtime_error("[Check] The false positives are not subsets of the extended syncmer difference");
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
                std::set_symmetric_difference(kmer_vectors[0].cbegin(), kmer_vectors[0].cend(), 
                                      kmer_vectors[1].cbegin(), kmer_vectors[1].cend(), 
                                      std::back_inserter(buffer));
            }
            emem::external_memory_vector<kmer_t> exact_sym_diff(max_ram, tmp_dir, "exact_diff");
            for (auto itr = exact_sym_diff.cbegin(); itr != exact_sym_diff.cend(); ++itr) {
                auto kmer = *itr;
                uint8_t* vptr = reinterpret_cast<uint8_t*>(&kmer);
                kmp::little2big(vptr, sizeof(kmer));
                exact_sym_diff.push_back(*reinterpret_cast<kmer_t*>(vptr));
            }
            if (exact_sym_diff.size() != true_sym_diff.size()) throw std::runtime_error("[Check] sizes do not match");
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
    parser.add_argument("extended-minuend")
        .help("extended IBLT from which to subtract");
    parser.add_argument("correction-minuend")
        .help("correction IBLT of the minuend (containing all k-mers)");
    parser.add_argument("subtrahend")
        .help("extended IBLT to subtract");
    parser.add_argument("correction-subtrahend")
        .help("correction IBLT of the subtrahend (containing all k-mers)");
    parser.add_argument("-l", "--list")
        .help("List IBLT to file. Use '.' for stdout");
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
    return parser;
}

} // namespace kmp