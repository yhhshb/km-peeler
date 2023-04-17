#include <cstdio>
#include <random>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <set>

#include <argparse/argparse.hpp>
#include "../../bundled/biolib/bundled/MurmurHash3.hpp"
#include "kmc_api/kmc_file.h"


struct hash_t {
	uint64_t msb;
	uint64_t lsb;
};

bool operator<(const hash_t& x, const hash_t& y) {
	return std::tie(x.msb, x.lsb) < std::tie(y.msb, y.lsb);
}

bool operator==(const hash_t& x, const hash_t& y) {
	return std::tie(x.msb, x.lsb) == std::tie(y.msb, y.lsb);
}

argparse::ArgumentParser get_wmh_parser()
{
	argparse::ArgumentParser parser("wmh");
	parser.add_description("Sketch a KMC database using WeightedMinHash (treating repeated elements as distinct values)");
	parser.add_argument("-i", "--input")
        .help("input kmc database")
		.required();
    parser.add_argument("-o", "--output")
        .help("output sketch")
		.required();
	parser.add_argument("-s", "--size")
		.help("number of elements to sample")
		.scan<'u', std::size_t>()
		.default_value(std::size_t(1000));
	parser.add_argument("--seed")
		.help("random seed")
		.scan<'u', uint32>()
		.default_value(uint32_t(42));
	return parser;
}

int wmh_main(const argparse::ArgumentParser& args)
{
	auto save_wmh_sketch = [](uint64_t seed, uint8_t k, std::vector<hash_t>& sketch, std::ofstream & ostrm)
	{
		ostrm.write(reinterpret_cast<char*>(&seed), sizeof(uint64_t));
		uint8_t hashdim = sizeof(hash_t);
		ostrm.write(reinterpret_cast<char*>(&hashdim), sizeof(uint8_t));
		ostrm.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
		seed = sketch.size();
		ostrm.write(reinterpret_cast<char*>(&seed), sizeof(uint64_t));
		ostrm.write(reinterpret_cast<char*>(sketch.data()), sketch.size() * sizeof(hash_t));
	};

	auto kmc_path = args.get<std::string>("--input");
	auto output_path = args.get<std::string>("--output");
	auto s = args.get<std::size_t>("--size");
	auto dummy_seed = args.get<uint32>("--seed");

	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_path)) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	CKMCFileInfo infos;
	kmcdb.Info(infos);
	unsigned int _kmer_length = infos.kmer_length;

	CKmerAPI kmer(_kmer_length);
	char skmer[_kmer_length];	
	uint32_t counter;
	hash_t hash;
	std::vector<hash_t> sketch;
	{
	std::set<hash_t> ssk;
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(skmer);
		for(std::uint32_t i = 0; i < counter; ++i)
		{
			MurmurHash3_x86_128(reinterpret_cast<void*>(skmer), _kmer_length, i, reinterpret_cast<void*>(&hash)); //get the seed for the column
			ssk.insert(hash);
			if(ssk.size() > s) ssk.erase(--(ssk.end()));
		}
	}
	kmcdb.Close();
	std::copy(ssk.cbegin(), ssk.cend(), std::back_inserter(sketch));
	}
	std::ofstream skstrm(output_path, std::ios_base::binary);
	save_wmh_sketch(dummy_seed, static_cast<uint8_t>(_kmer_length), sketch, skstrm);
	skstrm.close();
	return 0;
}

argparse::ArgumentParser get_compare_parser()
{
	argparse::ArgumentParser parser("compare");
	parser.add_description("Compare two WeightedMinHash sketches.\n\tOutput on stdout as a single float value (without newlines)");
	parser.add_argument("sketch1")
        .help("first sketch");
    parser.add_argument("sketch2")
        .help("second sketch");
	return parser;
}

int compare_main(const argparse::ArgumentParser& args)
{
	auto load_wmh_sketch = [](std::ifstream & istrm, uint64_t& seed, uint8_t& hashdim, uint8_t& k)
	{
		istrm.read(reinterpret_cast<char*>(&seed), sizeof(uint64_t));
		istrm.read(reinterpret_cast<char*>(&hashdim), sizeof(uint8_t));
		istrm.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
		uint64_t buffer;
		istrm.read(reinterpret_cast<char*>(&buffer), sizeof(decltype(buffer)));
		std::vector<hash_t> sketch(buffer);
		istrm.read(reinterpret_cast<char*>(sketch.data()), buffer * sizeof(hash_t));
		return sketch;
	};

	auto sketch1 = args.get<std::string>("sketch1");
	auto sketch2 = args.get<std::string>("sketch2");

	std::ifstream skstrm(sketch1, std::ios_base::binary);
	uint64_t seed1, seed2;
	uint8_t hashdim1, k1, hashdim2, k2;
	seed1 = hashdim1 = k1 = 1;
	seed2 = hashdim2 = k2 = 2;//just to be sure that they are read from the stream
	auto sk1 = load_wmh_sketch(skstrm, seed1, hashdim1, k1);
	skstrm.close();
	skstrm.open(sketch2, std::ios_base::binary);
	auto sk2 = load_wmh_sketch(skstrm, seed2, hashdim2, k2);

	if(seed1 != seed2 or hashdim1 != hashdim2 or k1 != k2) throw std::runtime_error("Incompatible sketches");

	std::size_t intersection = 0;
	std::size_t unon = 0;
	std::size_t seen = 0;
	auto hash_it1 = sk1.cbegin();
	auto hash_it2 = sk2.cbegin();
	while(seen < std::min(sk1.size(), sk2.size()) and (hash_it1 != sk1.cend() or hash_it2 != sk2.cend()))
	{
		if(hash_it1 != sk1.cend() and hash_it2 != sk2.cend() and *hash_it1 == *hash_it2)
		{
			unon += 1;
			intersection += 1;
			seen += 2;
			++hash_it1;
			++hash_it2;
		}
		else if (hash_it1 != sk1.cend() and (hash_it2 == sk2.cend() or *hash_it1 < *hash_it2))
		{
			unon += 1;
			seen += 1;
			++hash_it1;
		}
		else if (hash_it2 != sk2.cend() and (hash_it1 == sk1.cend() or *hash_it2 < *hash_it1))
		{
			unon += 1;
			seen += 1;
			++hash_it2;
		}
	}

	std::cout << static_cast<float>(static_cast<double>(intersection) / unon);
	return 0;
}

argparse::ArgumentParser get_i2cws_parser()
{
	argparse::ArgumentParser parser("i2cws");
	parser.add_description("Sketch a KMC database with I2CWS");
	parser.add_argument("-i", "--input")
        .help("input kmc database")
		.required();
    parser.add_argument("-o", "--output")
        .help("output sketch")
		.required();
	parser.add_argument("-s", "--size")
		.help("number of elements to sample")
		.scan<'u', std::size_t>()
		.default_value(std::size_t(1000));
	parser.add_argument("--seed")
		.help("random seed")
		.scan<'u', uint32>()
		.default_value(uint32_t(42));
	return parser;
}

int i2cws_main(const argparse::ArgumentParser& args)
{
	struct sample_t {
		double kstar;
		double ykstar;
	};

	auto save_i2cws_sketch = [](std::size_t seed, std::vector<sample_t>& sketch, std::ofstream & ostrm)
	{
		ostrm.write(reinterpret_cast<char*>(&seed), sizeof(decltype(seed)));
		seed = sketch.size();
		ostrm.write(reinterpret_cast<char*>(&seed), sizeof(decltype(seed)));
		ostrm.write(reinterpret_cast<char*>(sketch.data()), sketch.size() * sizeof(sample_t));
	};

	std::string kmc_path = args.get<std::string>("--input"); 
	auto output_path = args.get<std::string>("--output");
	auto s = args.get<std::size_t>("--size");
	auto rseed = args.get<uint32>("seed");
	
	std::vector<uint32_t> row_seeds(s);
	{
		std::mt19937_64 mersenne;
		std::uniform_int_distribution<uint32_t> sketch_dist(0, std::numeric_limits<uint32_t>::max());
		for(auto& seed : row_seeds) seed = sketch_dist(mersenne);
	}
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_path)) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	CKMCFileInfo infos;
	kmcdb.Info(infos);
	unsigned int _kmer_length = infos.kmer_length;

	CKmerAPI kmer(_kmer_length);
	char skmer[_kmer_length];	
	uint32_t counter;

	hash_t icws_seed; 
	sample_t dummy = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
	std::vector<sample_t> sketch(s, dummy);
	std::mt19937_64 gen;
	std::gamma_distribution<double> gamma(2, 1);
	std::uniform_real_distribution<double> uniform(0, 1);
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(skmer);
		for(std::size_t i = 0; i < s; ++i)
		{
			MurmurHash3_x64_128(reinterpret_cast<void*>(skmer), _kmer_length, row_seeds[i], reinterpret_cast<void*>(&icws_seed)); //get the seed for the column
			
			gen.seed(icws_seed.msb);
			double gamma1 = gamma(gen);
			double gamma2 = gamma(gen);
			double gamma3 = gamma(gen);
			gen.seed(icws_seed.lsb);
			double uniform1 = uniform(gen);
			double uniform2 = uniform(gen);

			double t_k2 = std::floor(std::log(counter) / gamma2 + uniform1);
			double z_k  = std::exp(gamma2 * (t_k2 - uniform2 + 1));
			double a_k  = gamma3 / z_k;
			
			if(a_k < sketch[i].kstar)
			{
				sketch[i].kstar = a_k;
				std::size_t t_kstar = std::floor(std::log(counter) / gamma1 + uniform1);
				sketch[i].ykstar = std::exp(gamma1 * (t_kstar - uniform1));
			}
		}
	}
	kmcdb.Close();
	
	std::ofstream skstrm(output_path, std::ios_base::binary);
	save_i2cws_sketch(rseed, sketch, skstrm);
	skstrm.close();
	return 0;
}

argparse::ArgumentParser get_full_parser()
{
	argparse::ArgumentParser parser("full");
	parser.add_description("Compute full Jaccard and full Weighted Jaccard between two KMC databases");
	parser.add_argument("kmc1")
        .help("first kmc database");
    parser.add_argument("kmc2")
        .help("second kmc database");
	return parser;
}

int full_main(const argparse::ArgumentParser& args)
{
	CKMCFile file1, file2;
	if (!file1.OpenForListing(args.get<std::string>("kmc1"))) throw std::runtime_error("Unable to open the first kmc database");
	if (!file2.OpenForListing(args.get<std::string>("kmc2"))) throw std::runtime_error("Unable to open the second kmc database");
	
	CKMCFileInfo infos;
	file1.Info(infos);
	unsigned int _kmer_length1 = infos.kmer_length;
	file2.Info(infos);
	unsigned int _kmer_length2 = infos.kmer_length;
	if (_kmer_length1 != _kmer_length2) {
		throw std::runtime_error("Error, the two data bases use different k-mer lengths");
	}

	CKmerAPI kmer1(_kmer_length1), kmer2(_kmer_length2);
	uint32_t counter1, counter2;
	counter1 = counter2 = 0;
	std::size_t numerator, denominator, unione, intersection;
	intersection = unione = numerator = denominator = 0;
	bool s2r1 = true;
	bool s2r2 = true;
	while(s2r1 || s2r2)
	{
		if(s2r1 && s2r2 && kmer1 == kmer2)
		{
			unione += 1;
			intersection += 1;
			numerator += std::min(counter1, counter2);
			denominator += std::max(counter1,counter2);
			s2r1 = file1.ReadNextKmer(kmer1, counter1);
			s2r2 = file2.ReadNextKmer(kmer2, counter2);
		}
		else if (s2r1 && (!s2r2 || kmer1 < kmer2))
		{
			unione += 1;
			denominator += counter1;
			s2r1 = file1.ReadNextKmer(kmer1, counter1);
		}
		else if (s2r2 && (!s2r1 || kmer2 < kmer1))
		{
			unione += 1;
			denominator += counter2;
			s2r2 = file2.ReadNextKmer(kmer2, counter2);
		}
		//std::cerr << "partial num = " << numerator << " partial den = " << denominator << std::endl;
	}
	if(intersection == 0 && unione == 0) intersection = unione = 1;
	if(numerator == 0 && denominator == 0) numerator = denominator = 1;
	std::cout << static_cast<double>(intersection)/static_cast<double>(unione) << 
	"," << 
	static_cast<double>(numerator)/static_cast<double>(denominator) <<
	"," << 
	unione - intersection;
	return 0;
}

argparse::ArgumentParser get_total_parser()
{
	argparse::ArgumentParser parser("total");
	parser.add_description("Get total number of k-mers in kmc database");
	parser.add_argument("kmcd")
        .help("a kmc database");
	return parser;
}

int total_main(const argparse::ArgumentParser& args) 
{
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(args.get<std::string>("kmcd"))) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	CKMCFileInfo infos;
	kmcdb.Info(infos);
	unsigned int _kmer_length = infos.kmer_length;

	CKmerAPI kmer(_kmer_length);	
	uint32_t counter;
	unsigned long long total = 0;

	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		total += counter;
	}

	std::cout << total;
	return 0;
}

int main(int argc, char* argv[])
{
	auto i2cws_parser = get_i2cws_parser();
	auto wmh_parser = get_wmh_parser();
	auto compare_parser = get_compare_parser();
	auto full_parser = get_full_parser();
	auto total_parser = get_total_parser();

	argparse::ArgumentParser program(argv[0]);
	program.add_subparser(i2cws_parser);
	program.add_subparser(wmh_parser);
	program.add_subparser(compare_parser);
	program.add_subparser(full_parser);
	program.add_subparser(total_parser);

	try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }

	if (program.is_subcommand_used(i2cws_parser)) return i2cws_main(i2cws_parser);
    else if (program.is_subcommand_used(wmh_parser)) return wmh_main(wmh_parser);
    else if (program.is_subcommand_used(compare_parser)) return compare_main(compare_parser);
	else if (program.is_subcommand_used(full_parser)) return full_main(full_parser);
	else if (program.is_subcommand_used(total_parser)) return total_main(total_parser);
    else std::cerr << program << std::endl;

	return 0;
}
