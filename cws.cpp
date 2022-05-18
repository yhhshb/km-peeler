#include <cstdio>
#include <random>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <set>
#include "./kmc_api/kmc_file.h"
extern "C" {
	#include "murmur3.h"
	#include "ketopt.h"
}

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

void print_subcommands()
{
	std::cerr << 
		"sketch\tsketch a KMC database with I2CWS\n" <<
		"compare\tcompare two I2CWS sketches\n" <<
		"full\tcompute full Jaccard and full Weighted Jaccard between two KMC databases\n" <<
		"total\tget total number of k-mers in kmc database\n" << 
		"\n" <<
		"To see the documentation for a given subcommand use it with option -h\n";
}

void print_full_help()
{
	std::cerr <<
		"the 'full' subcommand has two positional options only:\n"
		"[file1] first KMC database\n" <<
		"[file2] second KMC database\n" <<
		"\n" <<
		"output on stdout in the form <J>|<WJ> Example:\n" << 
		"0.3|0.4 -> 0.3 is the simple Jaccard while 0.4 is its weighted counterpart\n";
}

void print_sketch_help()
{
	std::cerr <<
		"sketch subcommand options:\n" <<
		"i\tinput KMC database to sketch (without extensions)\n" <<
		"o\toutput sketch file\n" <<
		"s\tnumber of elements in the sketch [1000]\n" <<
		"r\trandom seed [42]\n";
}

void print_compare_help()
{
	std::cerr <<
		"compare subcommand options:\n" <<
		"[sketch1] first sketch file\n" <<
		"[sketch2] second sketch file\n" <<
		"\n" <<
		"output on stdout in the form of <EWJ> Example:\n" <<
		"0.3 -> a single value (without newline or other characters)\n";
}

void print_total_help()
{
	std::cerr <<
		"the 'total' subcommand takes only one positional argument:\n"
		"[file1] a KMC database\n" <<
		"\n" <<
		"output on stdout the total number of k-mers in the database with repetitions (L1 norm)\n";
}

int wmh_main(int argc, char* argv[])
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

	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};
	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_path, output_path;
	std::size_t s = 1000;
	std::size_t dummy_seed = 42;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:h", longopts)) > 0)
	{
		if (c == 'r') {
			dummy_seed = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'i') {
			kmc_path = opt.arg;
		} else if (c == 's') {
			s = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'o') {
			output_path = opt.arg;
		} else {
			std::cerr << "Option (" << c << ") not available\n";
			print_sketch_help();
			return EXIT_FAILURE;
		}
	}
	if(kmc_path.length() == 0 or output_path.length() == 0) {
		if(kmc_path.length() == 0) std::cerr << "Option -i is mandatory\n";
		if(output_path.length() == 0) std::cerr << "Option -o is mandatory\n";
		print_sketch_help();
		return EXIT_FAILURE;
	}

	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_path)) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	// unsigned int _kmer_length;
	// unsigned int _mode;
	// unsigned int _counter_size;
	// unsigned int _lut_prefix_length;
	// unsigned int _signature_len;
	// unsigned int _min_count;
	// unsigned long long _max_count;
	// unsigned long long _total_kmers;
	// kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
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

int compare_main(int argc, char* argv[])
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

	if(argc != 3) 
	{
		print_compare_help();
		return EXIT_FAILURE;
	}

	std::ifstream skstrm(argv[1], std::ios_base::binary);
	uint64_t seed1, seed2;
	uint8_t hashdim1, k1, hashdim2, k2;
	seed1 = hashdim1 = k1 = 1;
	seed2 = hashdim2 = k2 = 2;//just to be sure that they are read from the stream
	auto sk1 = load_wmh_sketch(skstrm, seed1, hashdim1, k1);
	skstrm.close();
	skstrm.open(argv[2], std::ios_base::binary);
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

int i2cws_main(int argc, char* argv[])
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

	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};
	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_path, output_path;
	std::size_t s = 1000;
	std::size_t rseed = 42;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:h", longopts)) > 0)
	{
		if (c == 'r') {
			rseed = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'i') {
			kmc_path = opt.arg;
		} else if (c == 's') {
			s = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'o') {
			output_path = opt.arg;
		} else {
			std::cerr << "Option (" << c << ") not available\n";
			print_sketch_help();
			return EXIT_FAILURE;
		}
	}
	if(kmc_path.length() == 0 or output_path.length() == 0) {
		if(kmc_path.length() == 0) std::cerr << "Option -i is mandatory\n";
		if(output_path.length() == 0) std::cerr << "Option -o is mandatory\n";
		print_sketch_help();
		return EXIT_FAILURE;
	}
	
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
	
	// unsigned int _kmer_length;
	// unsigned int _mode;
	// unsigned int _counter_size;
	// unsigned int _lut_prefix_length;
	// unsigned int _signature_len;
	// unsigned int _min_count;
	// unsigned long long _max_count;
	// unsigned long long _total_kmers;
	// kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
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

int full_main(int argc, char* argv[])
{
	if(argc != 3) 
	{
		print_full_help();
		return EXIT_FAILURE;
	}

	CKMCFile file1, file2;
	if (!file1.OpenForListing(argv[1])) throw std::runtime_error("Unable to open the first kmc database");
	if (!file2.OpenForListing(argv[2])) throw std::runtime_error("Unable to open the second kmc database");
	
	// unsigned int _kmer_length1, _kmer_length2;
	// unsigned int _mode1, _mode2;
	// unsigned int _counter_size1, _counter_size2;
	// unsigned int _lut_prefix_length1, _lut_prefix_length2;
	// unsigned int _signature_len1, _signature_len2;
	// unsigned int _min_count1, _min_count2;
	// unsigned long long _max_count1, _max_count2;
	// unsigned long long _total_kmers1, _total_kmers2;
	// file1.Info(_kmer_length1, _mode1, _counter_size1, _lut_prefix_length1, _signature_len1, _min_count1, _max_count1, _total_kmers1);
	// file2.Info(_kmer_length2, _mode2, _counter_size2, _lut_prefix_length2, _signature_len2, _min_count2, _max_count2, _total_kmers2);
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

int get_total_kmers_main(int argc, char* argv[]) 
{
	if(argc != 2) 
	{
		print_total_help();
		return EXIT_FAILURE;
	}

	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(argv[1])) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	// unsigned int _kmer_length;
	// unsigned int _mode;
	// unsigned int _counter_size;
	// unsigned int _lut_prefix_length;
	// unsigned int _signature_len;
	// unsigned int _min_count;
	// unsigned long long _max_count;
	// unsigned long long _total_kmers;
	// kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
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
	ketopt_t om = KETOPT_INIT;
	int c;
	while((c = ketopt(&om, argc, argv, 0, "x", 0)) >= 0) {}
	if (om.ind == argc) {
		std::cerr << "[Error] Subcommand unavailable\n\n";
		print_subcommands();
		return 1;
	}

	if (std::strcmp(argv[om.ind], "sketch") == 0) {
		return wmh_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "compare") == 0) {
		return compare_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "full") == 0) {
		return full_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "total") == 0) {
		return get_total_kmers_main(argc - om.ind, &argv[om.ind]);
	} else {
		std::cerr << "Missing subcommand\n\n";
		print_subcommands();
	}
	return 0;
}
