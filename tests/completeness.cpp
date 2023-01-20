/**
 * Test completeness of retrieved symmetric difference using IBLTs by comparing to the exact result
 */

#include <zlib.h>
extern "C" {
    #include "../include/kseq.h"
}

#include <string>
#include <argparse/argparse.hpp>
#include "../include/kmer_view.hpp"

KSEQ_INIT(gzFile, gzread)

typedef uint64_t kmer_t;

int main(int argc, char* argv[])
{
    gzFile fp;
    kseq_t* seq;
    std::string input_filename;
    uint8_t k, m;

    argparse::ArgumentParser parser(argv[0]);
    parser.add_argument("first_fasta");
    parser.add_argument("second_fasta");
    parser.add_argument("k");
    parser.add_argument("z");
    parser.parse_args(argc, argv);
    
    auto input_filename = parser.get<std::string>("input_filename");

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) throw std::runtime_error("Unable to open the input file " + input_filename + "\n");
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        wrapper::kmer_view<kmer_t> view(seq->seq.s, seq->seq.l, k);
        for (std::size_t i = 0; i < seq->seq.l; ++i) {
            std::string dumb_kmer(&seq->seq.s[i], k);
            // TODO
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
}