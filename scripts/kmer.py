import builtins
import sys
import gzip
import math
import fastx

def smart_numeric_cast(s):
    def is_number(s: str):
        try:
            float(s)
            return True
        except ValueError:
            return False
    if is_number(s):
        n = float(s)
        if n.is_integer(): return int(s)#with large numbers casting float gives an approximation error, better to use the original string
        else: return n
    else:
        return s

def nibble2hex(msn: int, lsn: int):
    if msn >= 4 or lsn >= 4: raise ValueError("Bases must be 2bit encoded")
    n = msn * 4 + lsn
    return "{0:01X}".format(n)

def hex(km: str, enc_table: list) -> str:
    """Encode kmer as a 2bit-packed hex string"""
    en = [enc_table[ord(c)] for c in km]
    if len(en) % 2 == 1: en.append(0)
    l = len(en)
    res = list()
    print(en)
    for i in range(0, int(l / 2)): res.append(nibble2hex(en[i*2], en[i*2+1]))
    return "".join(res)

def canonical(kmer: str):
    kmer_rev = kmer.translate(fastx.comp_tab)[::-1]
    if kmer > kmer_rev: return kmer_rev
    return kmer

def set(seq: str, k: int, canon: bool, table: set[str]) -> set[str]:
    l = len(seq)
    if l < k: return
    if table == None: table = builtins.set()
    for i in range(l - k + 1):
        kmer = seq[i:(i+k)]
        ok = True
        for c in kmer: 
            if not (c in fastx.base_for): ok = False
        if ok: table.add(canonical(kmer) if canon else kmer)
    return table

def count(seq: str, k: int, canon: bool, table: dict) -> dict:
    l = len(seq)
    if l < k: return
    if table == None: table = dict()
    for i in range(l - k + 1):
        kmer = seq[i:(i+k)]
        ok = True
        for c in kmer: 
            if not (c in fastx.base_for): ok = False
        if ok:
            if canon: kmer = canonical(kmer)
            if kmer in table: table[kmer] += 1
            else: table[kmer] = 1
    return table

def diff(table1, table2, symmetric: bool = False) -> list:
    s1 = isinstance(table1, builtins.set)
    s2 = isinstance(table2, builtins.set)
    d1 = isinstance(table1, dict)
    d2 = isinstance(table2, dict)
    assert s1 or d1
    assert s2 or d2
    res = list()
    if s1: var1 = dict.fromkeys(table1, 1)
    else: var1 = table1
    if s2: var2 = dict.fromkeys(table2, 1)
    else: var2 = table2
    for key, val in var1.items():
        if key not in var2: 
            res.append(('i', key, val))
        else: 
            delta = val - var2[key]
            if delta: res.append(('i', key, delta))
    if symmetric:
        for key, val in var2.items():
            if key not in var1: 
                res.append(('j', key, val))
            else: 
                delta = val - var1[key]
                if delta: res.append(('j', key, delta))
    return res

def length(obj, sep=' ') -> int:
    import io
    if isinstance(obj, str):
        return len(obj)
    elif isinstance(obj, io.TextIOWrapper):
        kmer, _ = obj.readline().split(sep)
        obj.seek(0)
        return len(kmer)
    else:
        raise TypeError("Cannot handle object of type {}".format(type(obj)))

def minimizer(m: int, seq: str, hasher) -> str:
    l = len(seq)
    if l < m: return None
    minh = math.inf
    mini = None
    for i in range(l - m + 1):
        hval = hasher(seq[i:(i+m)].encode())
        if minh > hval: 
            minh = hval
            mini = i
    return seq[mini:(mini+m)], seq[:mini]+'|'+seq[mini+m:]

def check(seq: str, nt_table: list, bound: int):
    for c in seq:
        if nt_table[ord(c)] >= bound: return False
    return True

def jaccard(s1: builtins.set[str], s2: builtins.set[str]):
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    return intersection, union

class Spectrum:
    
    def __init__(self) -> "Spectrum":
        self.histogram = dict()
        
    def add(self, count: int):
        if count in self.histogram: self.histogram[count] += 1
        else: self.histogram[count] = 1

    def addFromDict(self, table: dict):
        for _, count in table.items():
            if count in self.histogram: self.histogram[count] += 1
            else: self.histogram[count] = 1

    def addFromFile(self, file, sep: str):
        if isinstance(file, str):
            file_handle = open(file, "r")
            close = True
        else: 
            file_handle = file
            close = False
        for line in file_handle:
            _, count = line.split(sep)
            count = smart_numeric_cast(count)
            if count in self.histogram: self.histogram[count] += 1
            else: self.histogram[count] = 1
        if close: file_handle.close()

    def empty(self) -> bool:
        return not bool(self.histogram)

    def L0Norm(self) -> int:
        return sum(list(self.histogram.values()))

    def L1Norm(self) -> int:
        L1 = 0
        for k, v in self.histogram.items(): L1 += smart_numeric_cast(k) * smart_numeric_cast(v)
        return L1

    def entropy(self) -> float:
        #L1 = self.L1Norm()
        L0 = self.L0Norm()
        H0 = 0
        for _, column in self.histogram.items():
            p = column/L0
            H0 += -p*math.log(p,2)
        return H0

    def getMaxColumn(self) -> int:
        return max(list(self.histogram.values()))

    def getMaxCount(self) -> int:
        mcounts = list()
        colmax = self.getMaxColumn()
        for count, column in self.histogram.items():
            if column == colmax: mcounts.append(count)
        mcounts.sort()
        return mcounts[0]

    def getOptimalEpsilon(self) -> float:
        """Dimension a Bloom Filter for the cold items of the spectrum.

        :return: The best epsilon for a Bloom Filter storing cold items only in order to minimize the false positives coming from the most common item.
        """
        L0 = self.L0Norm()
        N = self.getMaxColumn()
        return (L0 - N) / N
        #if N == None:
        #    self.getMaxCount()
        #    N = self._L1light
        #L1 = self.L1Norm()
        #return N / (L1 - N)

    def getOptimizedEpsilon(self, c_csf: float) -> float:
        L0 = self.L0Norm()
        N = self.getMaxColumn()
        c_bf = 1/math.log(2)
        epsilon = c_bf/c_csf * ((L0 - N)/N) * math.log(math.e, 2)
        if epsilon == 0: epsilon = math.inf #raise RuntimeError("epsilon = 0 | L0 = {} | N = {} | c_csf = {}".format(L0, N, c_csf))
        return epsilon

    def removeCount(self, count: int):
        if count in self.histogram: 
            del self.histogram[count]

def set_main(args):
    table = builtins.set()
    for f in args.i:
        if f.endswith(".gz"): fi = gzip.open(f, "rt")
        else: fi = open(f, "r")
        for _, seq, _ in fastx.read(fi):
            set(seq, args.k, args.c, table)
        fi.close()
    if not args.i:
        for _, seq, _ in fastx.read(sys.stdin):
            set(seq, args.k, args.c, table)
    if (args.o): fo = open(args.o, "w")
    else: fo = sys.stdout
    for k in table: fo.write("{}\n".format(k))
    fo.close()

def count_main(args):
    assert len(args.sep) == 1
    table = dict()
    for f in args.i:
        if f.endswith(".gz"): fi = gzip.open(f, "rt")
        else: fi = open(f, "r")
        for _, seq, _ in fastx.read(fi):
            count(seq, args.k, args.c, table)
        fi.close()
    if not args.i:
        for _, seq, _ in fastx.read(sys.stdin):
            count(seq, args.k, args.c, table)
    if (args.o): fo = open(args.o, "w")
    else: fo = sys.stdout
    for k, v in table.items(): fo.write("{}{}{}\n".format(k, args.sep, v))
    fo.close()

def diff_main(args):
    table1 = dict()
    table2 = dict()
    with open(args.i, "r") as th:
        for _, seq, _ in fastx.read(th):
            count(seq, args.k, args.c, table1)
    with open(args.j, "r") as th:
        for _, seq, _ in fastx.read(th):
            count(seq, args.k, args.c, table2)
    difference = diff(table1, table2, args.s)
    for s, km, delta in difference:
        sys.stdout.write("{},{},{}\n".format(s, km, delta))

def tab_diff_main(args):
    assert len(args.sep) == 1
    si = builtins.set()
    with open(args.i, "r") as ih:
        for line in ih:
            km = line.strip().split(args.sep)[0]
            if args.g:
                if check(km, fastx.enc_nuc_table, 4): si.add(km)
            else:
                si.add(km)
    sj = builtins.set()
    with open(args.j, "r") as jh:
        for line in jh:
            km = line.strip().split(args.sep)[0]
            if args.g:
                if check(km, fastx.enc_nuc_table, 4): sj.add(km)
            else:
                sj.add(km)
    difference = diff(si, sj, args.s)
    for s, km, delta in difference:
        sys.stdout.write("{},{},{}\n".format(s, km, delta))

def tab_jaccard_main(args):
    assert len(args.sep) == 1
    si = builtins.set()
    with open(args.i, "r") as ih:
        for line in ih:
            km = line.strip().split(args.sep)[0]
            if args.g: 
                if check(km, fastx.enc_nuc_table, 4): si.add(km)
            else: 
                si.add(km)
    sj = builtins.set()
    with open(args.j, "r") as jh:
        for line in jh:
            km = line.strip().split(args.sep)[0]
            if args.g: 
                if check(km, fastx.enc_nuc_table, 4): sj.add(km)
            else: 
                sj.add(km)
    i, u = jaccard(si, sj)
    sys.stderr.write("{}/{}\n".format(i, u))
    print(i / u)

def histogram_main(args):
    assert len(args.sep) == 1
    sp = Spectrum()
    if args.i: ih = open(args.i, "r")
    else: ih = sys.stdin
    for line in ih:
        sp.add(int(line.strip().split(args.sep)[1]))
    ih.close()
    if args.o: oh = open(args.o, "w")
    else: oh = sys.stdout
    for v, c in sorted([(v,c) for v,c in sp.histogram.items()]):
        oh.write("{},{}\n".format(v,c))
    oh.close()
    sys.stderr.write("entropy = {}\n".format(sp.entropy()))

def syncmers_main(args):
    import xxhash
    si = builtins.set()
    if args.i: ih = open(args.i, "r")
    else: ih = sys.stdin
    for line in ih:
        km = line.strip().split(args.sep)[0]
        if args.g:
            if check(km, fastx.enc_nuc_table, 4): si.add(km)
        else:
            si.add(km)
    ih.close()
    if args.o: fo = open(args.o, "w")
    else: fo = sys.stdout
    for km in si:
        rem = minimizer(args.z, km, xxhash.xxh32_intdigest)[1]
        if (rem.startswith('|') or rem.endswith('|')):
            fo.write("{}\n".format(km))
    fo.close()

def tofasta_main(args):
    counter = 0
    if args.i: ih = open(args.i, "r")
    else: ih = sys.stdin
    if args.o: oh = open(args.o, "w")
    else: oh = sys.stdout
    for kmer in ih:
        oh.write(">{}\n{}\n".format(counter, kmer.strip()))
        counter += 1
    oh.close()
    ih.close()

def main(args):
    if (args.command == "set"): set_main(args)
    elif (args.command == "count"): count_main(args)
    elif (args.command == "diff"): diff_main(args)
    elif (args.command == "tabdiff"): tab_diff_main(args)
    elif (args.command == "jaccard"): tab_jaccard_main(args)
    elif (args.command == "histogram"): histogram_main(args)
    elif (args.command == "syncmers"): syncmers_main(args)
    elif (args.command == "2fasta"): tofasta_main(args)
    else: parser.print_help(sys.stderr)
    
def setup_parser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_set = subparsers.add_parser("set", help="Compute set of k-mers")
    parser_set.add_argument("-i", help="input FASTX file [stdin]", type=str, nargs='+', default=[])
    parser_set.add_argument("-o", help="output (one k-mer per line) [stdout]", type=str)
    parser_set.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_set.add_argument("-c", help="canonical k-mers", action="store_true")

    parser_count = subparsers.add_parser("count", help="Count k-mers")
    parser_count.add_argument("-i", help="input files (fasta or fastq) [stdin]", type=str, nargs='+', default=[])
    parser_count.add_argument("-o", help="output count table [stdout]", type=str)
    parser_count.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_count.add_argument("-c", help="canonical k-mers", action="store_true")
    parser_count.add_argument("--sep", help="separator between kmers and their counts", type=str, default=',')

    parser_diff = subparsers.add_parser("diff", help="Compute k-mer difference between two fastx files")
    parser_diff.add_argument("-i", help="first fasta file", type=str, required=True)
    parser_diff.add_argument("-j", help="second fasta file", type=str, required=True)
    parser_diff.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_diff.add_argument("-c", help="canonical k-mers", action="store_true")
    parser_diff.add_argument("-s", help="symmetric difference", action="store_true")

    parser_tabdiff = subparsers.add_parser("tabdiff", help="Compute k-mer difference between two k-mer tables")
    parser_tabdiff.add_argument("-i", help="first table", type=str, required=True)
    parser_tabdiff.add_argument("-j", help="second table", type=str, required=True)
    parser_tabdiff.add_argument("-s", help="symmetric difference", action="store_true")
    parser_tabdiff.add_argument("-g", help="ignore strings with letters other than [ACGT]", action="store_true")
    parser_tabdiff.add_argument("--sep", help="separator between kmers and their counts", type=str, default=',')

    parser_jaccard = subparsers.add_parser("jaccard", help="Compute jaccard index between two k-mer tables")
    parser_jaccard.add_argument("-i", help="first table", type=str, required=True)
    parser_jaccard.add_argument("-j", help="second table", type=str, required=True)
    parser_jaccard.add_argument("-g", help="ignore strings with letters other than [ACGT]", action="store_true")
    parser_jaccard.add_argument("--sep", help="separator between kmers and their counts", type=str, default=',')

    parser_histogram = subparsers.add_parser("histogram", help="Compute histogram from k-mer count table and relative statistics")
    parser_histogram.add_argument("-i", help="k-mer count table [stdin]", type=str)
    parser_histogram.add_argument("-o", help="output file [stdout]", type=str)
    parser_histogram.add_argument("--sep", help="separator between kmers and their counts [,]", type=str, default=',')

    parser_syncmers = subparsers.add_parser("syncmers", help="Subsample a given input table by selecting syncmers")
    parser_syncmers.add_argument("-i", help="k-mer table [stdin]", type=str, required=False)
    parser_syncmers.add_argument("-o", help="output (one syncmer per line) [stdout]", type=str)
    parser_syncmers.add_argument("-z", help="syncmer parameter (must be less than the k-mer length)", type=int, required=True)
    parser_syncmers.add_argument("-g", help="ignore strings with letters other than [ACGT]", action="store_true")
    parser_syncmers.add_argument("--sep", help="separator between kmers and their counts", type=str, default=',')

    parser_2fasta = subparsers.add_parser("2fasta", help="transform a stream of k-mers into a fasta-formatted stream")
    parser_2fasta.add_argument("-i", help="k-mer table [stdin]", type=str, required=False)
    parser_2fasta.add_argument("-o", help="output (one syncmer per line) [stdout]", type=str, required=False)

    return parser

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args(sys.argv)
    main(args)

