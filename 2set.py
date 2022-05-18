import sys

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def main(args):
    if args.i: th = open(args.i, "r")
    else: th = sys.stdin
    unique = set()
    for line in th:
        line = line.strip()
        #if th: 
        if args.c:
            kmer_rev = line.translate(comp_tab)[::-1]
            if line > kmer_rev: line = kmer_rev
        unique.add(line)
    th.close()
    if args.o: oh = open(args.o, "w")
    else: oh = sys.stdout
    for s in unique:
        oh.write("{}\n".format(s))
    oh.close()

def setup_parser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input stream of fragments [stdin]", type=str)
    parser.add_argument("-o", help="output unique fragments [stdout]", type=str)
    parser.add_argument("-c", help="canonical k-mers only", action="store_true")

    return parser

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    main(args)