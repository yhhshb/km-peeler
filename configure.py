import os
import sys
import subprocess

compile_options_filename = "compile_options.h"

def configure_length(algorithm: str, k: int, m: int):
    assert k > 0
    assert 0 <= m <= k
    if algorithm == "kmers": return k
    elif algorithm == "minimizers": return 2*k-m #space for k-mers grouped together
    elif algorithm == "syncmers": return k #allocate space to store syncmers
    elif algorithm == "segmentation": return 2*k #allocate enough space for the segmentation algorithm
    #elif algorithm == "min-hash-collection": return int((k + 7) / 8) * m # = CEILING(hash_width, 8) * number_of_hashes_in_minhash_sketch
    else: raise ValueError("This should never happen because the parser should catch it with options")

def configure_for_fixed_length_sequences(path: str, algorithm: str, k: int, z: int, always_make: bool, debug: bool):
    assert 0 <= k
    assert 0 <= z < k
    aldiff_path = path #same folder as source code
    compile_options_path = os.path.join(aldiff_path, compile_options_filename)
    with open(compile_options_path, "w") as ch:
        ch.write("#ifndef COMPILE_OPTIONS_H\n")
        ch.write("#define COMPILE_OPTIONS_H\n")
        ch.write("\n")
        # ch.write("#define READ_KMERS         /*Read a set of k-mers*/\n")
        dnalen = configure_length(algorithm, k, z)
        print(algorithm)
        print("-->", dnalen)
        ch.write("#define STORE_SEQUENCES {}\n".format(dnalen))
        ch.write("\n")
        ch.write("{}#define DEBUG{}\n".format("" if debug else "/*", "" if debug else "*/"))
        ch.write("\n")
        ch.write("#endif\n")
    out = subprocess.run(["cd {} && make{}".format(aldiff_path, " --always-make" if always_make else "")], capture_output=False, shell=True) #--always-make
    if out.returncode != 0:
        sys.stderr.write("Unable to re-build aldiff for fixed-length sequences\n")
        exit(out.returncode)

def configure_for_hashes(path: str, bit_width: int, minhash_size: int, always_make: bool, debug: bool):
    assert 0 < bit_width <= 64
    aldiff_path = path #same folder as source code
    compile_options_path = os.path.join(aldiff_path, compile_options_filename)
    hash_width_in_bits = int((bit_width + 7) / 8) * minhash_size * 8
    with open(compile_options_path, "w") as ch:
        ch.write("#ifndef COMPILE_OPTIONS_H\n")
        ch.write("#define COMPILE_OPTIONS_H\n")
        ch.write("\n")
        # ch.write("#define READ_KMERS         /*Read a set of k-mers*/\n")
        ch.write("#define STORE_HASHES {}\n".format(hash_width_in_bits))
        ch.write("\n")
        ch.write("{}#define DEBUG{}\n".format("" if debug else "/*", "" if debug else "*/"))
        ch.write("\n")
        ch.write("#endif\n")
    out = subprocess.run(["cd {} && make{}".format(aldiff_path, " --always-make" if always_make else "")], capture_output=False, shell=True)
    if out.returncode != 0:
        sys.stderr.write("Unable to re-build aldiff for hashes\n")
        exit(out.returncode)

def configure_for_variable_length_sequences(path: str, maximum_length: int, always_make: bool, debug: bool):
    assert maximum_length >= 0
    aldiff_path = path
    compile_options_path = os.path.join(aldiff_path, compile_options_filename)
    with open(compile_options_path, "w") as ch:
        ch.write("#ifndef COMPILE_OPTIONS_H\n")
        ch.write("#define COMPILE_OPTIONS_H\n")
        ch.write("\n")
        # ch.write("#define READ_KMERS         /*Read a set of k-mers*/\n")
        ch.write("#define STORE_VLSEQUENCES {}\n".format(maximum_length))
        ch.write("\n")
        ch.write("{}#define DEBUG{}\n".format("" if debug else "/*", "" if debug else "*/"))
        ch.write("\n")
        ch.write("#endif\n")
    out = subprocess.run(["cd {} && make{}".format(aldiff_path, " --always-make" if always_make else "")], capture_output=False, shell=True)
    if out.returncode != 0:
        sys.stderr.write("Unable to re-build aldiff for variable length fragments\n")
        exit(out.returncode)

def main(args):
    if args.command == "sequences": return configure_for_fixed_length_sequences(args.__exepath, args.algorithm, args.k, args.z, args.always, args.debug)
    elif args.command == "hashes": return configure_for_hashes(args.__exepath, args.l, args.s, args.always, args.debug)
    elif args.command == "vls": return configure_for_variable_length_sequences(args.__exepath, args.l, args.always, args.debug)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_sequences = subparsers.add_parser("sequences", help="Configure aldiff to store sequences for the values of k and m used by its subcommand")
    parser_sequences.add_argument("-a", "--algorithm", help="algorithm to compute the set (minimizers, syncmers) [minimizers]", type=str, default="minimizers", choices=["minimizers", "syncmers", "segmentation"])
    parser_sequences.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_sequences.add_argument("-z", help="minimizer length used for improving space requirements (see algorithm for actual usage)", type=int, required=True)
    parser_sequences.add_argument("--debug", help="activate debugging information to stderr", action="store_true")
    parser_sequences.add_argument("--always", help="recompile everything", action="store_true")

    parser_hashes = subparsers.add_parser("hashes", help="Configure aldiff to store hashes")
    parser_hashes.add_argument("-l", help="hash length in bits [64]", type=int, default=64)
    parser_hashes.add_argument("-s", help="store minHash sketches (one sketch is like a very big hash value)", type=int, default=1)
    parser_hashes.add_argument("--debug", help="activate debugging information to stderr", action="store_true")
    parser_hashes.add_argument("--always", help="recompile everything", action="store_true")

    parser_vls = subparsers.add_parser("vls", help="Configure aldiff to store variable length sequences")
    parser_vls.add_argument("-l", help="Maximum length of the sequences to be stored inside an IBF", type=int, required=True)
    parser_vls.add_argument("--debug", help="activate debugging information to stderr", action="store_true")
    parser_vls.add_argument("--always", help="recompile everything", action="store_true")
    
    return parser

if __name__ == "__main__":
    mypath = os.path.dirname(sys.argv[0])
    abspath = os.path.abspath(mypath)
    parser = parser_init()
    args = parser.parse_args(sys.argv)
    setattr(args, "__exepath", abspath)
    main(args)