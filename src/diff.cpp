#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include <array>
#include "../include/diff.hpp"
#include "../include/compile_constants.hpp"
#include "../include/external_memory_vector.hpp"
#include "../include/kmer_view.hpp"
#include "../include/syncmer_sampler.hpp"
#include "../include/io.hpp"
#include "../include/IBLT.hpp"
#include "../include/utils.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kmp {

int diff_main(const argparse::ArgumentParser& args)
{
    std::size_t loaded_bytes, m_size, s_size;
    uint64_t max_ram = args.get<uint64_t>("--max-ram");
    std::string tmp_dir = args.get<std::string>("--tmp-dir");
    std::string m_filename = args.get<std::string>("minuend");
    std::string s_filename = args.get<std::string>("subtrahend");
    std::string jaccard_filename, list_filename;
    IBLT mned = load(m_filename, loaded_bytes);
    std::cerr << "[Info] minuend: loaded " << loaded_bytes << " bytes\n";
    m_size = mned.size();
    // mned.print_config(std::cerr);
    { // memory management of sbhd
        IBLT sbhd = load(s_filename, loaded_bytes);
        s_size = sbhd.size();
        std::cerr << "[Info] subtrahend: loaded " << loaded_bytes << " bytes\n";
        // sbhd.print_config(std::cerr);
        mned.subtract(sbhd);
    }
    if (args.is_used("--output")) io::store(mned, args.get<std::string>("--output"));
    if (args.is_used("--jaccard")) jaccard_filename = args.get<std::string>("--jaccard");
    if (args.is_used("--list")) list_filename = args.get<std::string>("--list");
    emem::external_memory_vector<kmer_t> symmetric_difference(max_ram, tmp_dir, "check");

    if (args.is_used("--check")) {
        gzFile fp;
        kseq_t* seq;
        uint8_t z, offset1, offset2;
        auto files = args.get<std::vector<std::string>>("--check");
        if (files.size() != 2) throw std::invalid_argument("check option requires two arguments");
        if (args.is_used("-z")) z = args.get<uint8_t>("-z");
        else throw std::invalid_argument("syncmer parameter z used to build the IBLTs is required when checking");
        if (args.is_used("--offset")) {
            offset1 = args.get<uint8_t>("--offset");
            assert(mned.get_k() > z);
            offset2 = offset1 ? offset1 : mned.get_k() - z;
        } else {
            throw std::invalid_argument("offset option used to build the IBLTs is required when checking");
        }
        std::array<emem::external_memory_vector<kmer_t>, 2> kmer_vectors = {
            emem::external_memory_vector<kmer_t> (max_ram, tmp_dir, "set1"),
            emem::external_memory_vector<kmer_t> (max_ram, tmp_dir, "set2")
        };
        hash::minimizer_position_extractor mmp_extractor(mned.get_k(), z);
        for (std::size_t i = 0; i < files.size(); ++i) {
            auto& f = files.at(i);
            if ((fp = gzopen(f.c_str(), "r")) == NULL) throw std::runtime_error("Unable to open the input file " + f + "\n");
            seq = kseq_init(fp);
            while (kseq_read(seq) >= 0) {
                wrapper::kmer_view<kmer_t> view(seq->seq.s, seq->seq.l, mned.get_k());
                sampler::syncmer_sampler syncmers(view.cbegin(), view.cend(), mmp_extractor, offset1, offset2);
                for (auto itr = syncmers.cbegin(); itr != syncmers.cend(); ++itr) {
                    kmer_vectors[i].push_back(*itr);
                }
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
        }
        std::cerr << "[1] checkpoint" << std::endl;
        std::set_symmetric_difference(kmer_vectors[0].cbegin(), kmer_vectors[0].cend(), 
                                      kmer_vectors[1].cbegin(), kmer_vectors[1].cend(), 
                                      std::back_inserter(symmetric_difference));
        std::cerr << "[2] checkpoint" << std::endl;
    }
    std::vector<std::vector<uint8_t>> p, n;
    IBLT::failure_t err;
    err = mned.list(p, n);
    if (list_filename.length()) {
        auto cout_buffer_save = std::cout.rdbuf();
        std::ofstream myostrm;
        if (list_filename != ".") { // redirect cout
            myostrm.open(list_filename);
            std::cout.rdbuf(myostrm.rdbuf());
        }

        if (not err) {
            for (auto& v : p) std::cout << vec2kmer(v, mned.get_k()) << "\n";
            for (auto& v : n) std::cout << vec2kmer(v, mned.get_k()) << "\n";
        } else {
            std::cout << err << "\n";
        }

        std::cout.rdbuf(cout_buffer_save); // restore cout
    }

    if (jaccard_filename.length()) {
        std::size_t m_diff_size = p.size(); 
        std::size_t s_diff_size = n.size();
        auto cout_buffer_save = std::cout.rdbuf();
        std::ofstream myostrm;
        if (jaccard_filename != ".") { // redirect cout
            myostrm.open(jaccard_filename);
            std::cout.rdbuf(myostrm.rdbuf());
        }
        std::cout << m_filename << "," << s_filename << ",";
        auto intersection_size = m_size - m_diff_size;
        assert(intersection_size == s_size - s_diff_size);
        double j = double(intersection_size) / double(m_size + s_size - intersection_size);
        if (not err) std::cout << j;
        else std::cout << err;
        std::cout << "\n";
        std::cout.rdbuf(cout_buffer_save); // restore cout
    }

    if (args.is_used("--check")) {
        if (not err) {
            std::vector<kmer_t> iblt_diff;
            iblt_diff.reserve(symmetric_difference.size());
            for (auto& v : p) {
                uint8_t* vptr = v.data();
                kmp::little2big(vptr, v.size());
                iblt_diff.push_back(*reinterpret_cast<kmer_t*>(vptr));
            }
            for (auto& v : n) {
                uint8_t* vptr = v.data();
                kmp::little2big(vptr, v.size());
                iblt_diff.push_back(*reinterpret_cast<kmer_t*>(vptr));
            }
            std::sort(iblt_diff.begin(), iblt_diff.end());
            std::vector<kmer_t> dummy;
            std::set_symmetric_difference(symmetric_difference.cbegin(), symmetric_difference.cend(),
                                            iblt_diff.cbegin(), iblt_diff.cend(), std::back_inserter(dummy));
            std::cerr << "Peelable sketch check: ";
            if (dummy.size()) std::cerr << "FAILED";
            else std::cerr << "OK";
            std::cerr << std::endl;
        } else {
            std::cerr << "Unable to check due to unpeelable sketch\n";
        }
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
    parser.add_argument("-c", "--check")
        .help("Original fasta files used to build the two IBLTs (only syncmers are supported)")
        .nargs(2);
    parser.add_argument("-z")
        .help("syncmer parameter")
        .scan<'u', uint8_t>();
    parser.add_argument("-f", "--offset")
        .help("z-mer offset inside syncmers")
        .scan<'u', uint8_t>()
        .default_value(uint8_t(0));
    parser.add_argument("-m", "--max-ram")
        .help("RAM limit")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(6));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

}