CFLAGS     = -pedantic -Wall -Werror -Wno-format -O3 -DNDEBUG=1 #-g
CXXFLAGS   = -std=c++14 -pedantic -Wall -Werror -Wno-char-subscripts -Wno-vla-extension -O3
# ASAN_FLAGS = -fsanitize=address -fno-omit-frame-pointer -Wno-format-security
# ASAN_LIBS  = -static-libasan

.PHONY:clean

all: ibltseq cws

ibltseq: aldiff.o kmers_main.o minimizers_main.o syncmers_main.o sample_main.o build_main.o diff_main.o list_main.o jaccard_main.o collection_main.o minHash.o print_main.o dump_main.o ibflib.o mmlib.o constants.o err.o endian_fixer.o kalloc.o murmur3.o
	$(CC) $(CFLAGS) $(ASAN_FLAGS) $(ASAN_LIBS) -o $@ $^ -lm -lz

aldiff.o: aldiff.c kmers_main.h minimizers_main.h syncmers_main.h sample_main.h build_main.h diff_main.h dump_main.h list_main.h print_main.h mmlib.h constants.h err.h kvec2.h kseq.h ketopt.h
	$(CC) $(CFLAGS) -c aldiff.c

kmers_main.o: kmers_main.h kmers_main.c constants.h err.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -c kmers_main.c

syncmers_main.o: syncmers_main.h syncmers_main.c err.h mmlib.o ketopt.h kvec2.h kseq.h
	$(CC) $(CFLAGS) -c syncmers_main.c

minimizers_main.o: minimizers_main.h minimizers_main.c err.h mmlib.o ketopt.h kvec2.h kseq.h
	$(CC) $(CFLAGS) -c minimizers_main.c

sample_main.o: sample_main.h sample_main.c err.h ketopt.h murmur3.h
	$(CC) $(CFLAGS) -c sample_main.c

build_main.o: build_main.h build_main.c err.h ibflib.o ketopt.h
	$(CC) $(CFLAGS) -c build_main.c

diff_main.o: diff_main.h diff_main.c ibflib.h err.h ketopt.h
	$(CC) $(CFLAGS) -c diff_main.c

list_main.o: list_main.h list_main.c ibflib.o err.h ketopt.h
	$(CC) $(CFLAGS) -c list_main.c

jaccard_main.o: jaccard_main.h jaccard_main.c ibflib.h err.h ketopt.h
	$(CC) $(CFLAGS) -c jaccard_main.c

collection_main.o: collection_main.h collection_main.c ibflib.h err.h constants.o
	$(CC) $(CFLAGS) -c collection_main.c

print_main.o: print_main.h print_main.c ibflib.h err.h ketopt.h
	$(CC) $(CFLAGS) -c print_main.c

dump_main.o: dump_main.h dump_main.c ibflib.h err.h ketopt.h
	$(CC) $(CFLAGS) -c dump_main.c

ibflib.o: ibflib.h ibflib.c endian_fixer.h err.h constants.o murmur3.h
	$(CC) $(CFLAGS) -c ibflib.c

mmlib.o: mmlib.h mmlib.c constants.o kalloc.o kvec2.h 
	$(CC) $(CFLAGS) -c mmlib.c

minHash.o: minHash.h minHash.c endian_fixer.h constants.h err.h murmur3.h
	$(CC) $(CFLAGS) -c minHash.c

constants.o: constants.h constants.c err.h compile_options.h
	$(CC) $(CFLAGS) -c constants.c

err.o: err.h err.c
	$(CC) $(CFLAGS) -c err.c

endian_fixer.o: endian_fixer.h endian_fixer.c
	$(CC) $(CFLAGS) -c endian_fixer.c

kalloc.o: kalloc.h kalloc.c
	$(CC) $(CFLAGS) -c kalloc.c

cws: cws.cpp murmur3.o
	$(CXX) $(CXXFLAGS) $(ASAN_FLAGS) $(ASAN_LIBS) -o $@ cws.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp murmur3.o

murmur3.o: murmur3.h murmur3.c
	$(CC) $(CFLAGS) -c murmur3.c

clean:
	rm -f *.o
	rm -f ibltseq
	rm -f cws
	