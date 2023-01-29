import os
import sys
import pathlib
import random
import kmc_wrapper
import mash_wrapper
import kmp_wrapper
import essentials
import pandas as pd

jaccard_choices = ["exact", "minHash", "iblt"]

def jaccard_experiment_main(args):
    '''Compute Jaccard using one of multiple techniques: {}

    Output is in CSV format: 
    - keys are <reference,query>
    - KMC: <>
    '''.format(jaccard_choices)
    kmp_executable = args.__bin_path.joinpath("kmp")
    cws_executable = args.__bin_path.joinpath("cws")
    assert kmp_executable.exists()
    assert 0 < args.k <= 32
    assert 0 <= args.z <= args.k
    assert 3 <= args.repetitions <= 7
    assert (args.query_number != None and args.query_number >= 0) or (args.query_fraction != None and 0 <= args.query_fraction <= 1)
    assert args.sketch_size >= 0
    assert args.input_folder.exists() and args.input_folder.is_dir()
    assert args.tmp_dir.exists() and args.tmp_dir.is_dir()
    wdir = args.tmp_dir

    _, ddir = essentials.get_path_leaf(args.input_folder)
    random.seed(args.seed)

    tmp_mash_file = wdir.joinpath(essentials.get_random_name("csv"))
    tmp_kmp_file = wdir.joinpath(essentials.get_random_name("csv"))
    tmp_kmc_file = wdir.joinpath(essentials.get_random_name("csv"))
    kmc_dir = wdir.joinpath("kmc")
    kmc_dataset_dir = kmc_dir.joinpath(ddir) # used later to retrieve kmc databases
    kmc_tmp_dir = wdir.joinpath("kmc_tmp")

    mash_dir = wdir.joinpath("mash")
    mash_ref = mash_dir.joinpath(essentials.get_random_name(""))

    kmp_dir = wdir.joinpath("kmp")
    extkmp_dir = wdir.joinpath("extended_syncmers_ibf")

    os.makedirs(kmc_dir, exist_ok=True)
    os.makedirs(kmc_tmp_dir, exist_ok=True)
    os.makedirs(mash_dir, exist_ok=True)
    os.makedirs(kmp_dir, exist_ok=True)
    os.makedirs(extkmp_dir, exist_ok=True)

    k = args.k
    z = args.z
    r = args.repetitions
    eps = args.epsilon
    s = mash_wrapper.get_size_from_byte_size(args.sketch_size, k) # args.sketch_size is in Bytes
    n = kmp_wrapper.get_maximum_difference_from_byte_size(args.sketch_size, k, r, eps)
    if s == 0: raise RuntimeError("minHash sketches of size 0 given the current parameters")
    if n == 0: raise RuntimeError("IBF sketches of size 0 given the current parameters")
    # sys.stderr.write("mash s = {}, kmp n = {} (--> sketch size = {})\n".format(s, n, kmp_wrapper.get_bit_size(n, k, r, eps)/8))

    # partition file set into query (Q) and reference (R)
    data_set = essentials.get_fastx_in_folder(args.input_folder)
    if args.query_number != None and args.query_fraction == None: query_size = args.query_number
    elif args.query_number == None and args.query_fraction != None: query_size = max(1, int(len(data_set) * args.query_fraction))
    else: raise RuntimeError("Something went wrong with the mutually-exclusive options 'qnumber' and 'qfraction'")
    if not query_size: query_size = len(data_set)
    query_set = sorted(random.sample(data_set, query_size))
    if (query_size == len(data_set)): reference_set = data_set
    else: reference_set = sorted(list(set(data_set) - set(query_set)))

    assert query_size

    # Count files using KMC
    sys.stderr.write("Build KMC databases and compute pairwise exact Jaccard\n")
    _ = kmc_wrapper.kmc_count_folder(args.input_folder, kmc_dir, [k], kmc_tmp_dir, args.max_ram, False, args.canonical, args.force_kmc)
    if args.force_kmc or not pathlib.Path(tmp_kmc_file).exists():
        with open(tmp_kmc_file, "w") as kmchandle:
            kmchandle.write("reference,query,exj\n")
            for query in query_set:
                for reference in reference_set:
                    _, _, _, _, refkmc = kmc_wrapper.get_kmc_paths(k, reference, kmc_dataset_dir) #retrieve reference file name
                    _, _, _, _, qrykmc = kmc_wrapper.get_kmc_paths(k, query, kmc_dataset_dir) #retrieve query file name
                    exj, _, _ = kmc_wrapper.kmc_jaccard(cws_executable, qrykmc, refkmc) # ignore weighted jaccard and true symmetric difference size
                    kmchandle.write("{},{},{}\n".format(reference, query, exj))

    # Build MASH sketches on input sequences
    sys.stderr.write("Mash estimations\n")
    mash_wrapper.mash_sketch(reference_set, mash_ref, k, s, args.seed, args.canonical, args.force_mash)
    mash_wrapper.mash_dist(mash_ref.with_suffix(".msh"), query_set, k, s, args.seed, tmp_mash_file, args.force_mash)

    # sketch using kmp
    sys.stderr.write("Build kmp sketches and compute Jaccard approximations\n")
    header_names = list()
    for sampling_rate in args.sampling_rates:
        assert sampling_rate >= 0
        header_names.append("syncj r={}".format(sampling_rate))
        kmp_wrapper.sketch_folder(kmp_executable, args.input_folder, kmp_dir, n, k, z, r, eps, args.seed, wdir, args.max_ram, args.force_kmp) #, args.canonical)
    
    # compute similarity scores between Q and R using [exact Jaccard, minHash, syncmers + IBF]
    if args.force_kmp or not pathlib.Path(tmp_kmp_file).exists():
        with open(tmp_kmp_file, "w") as kmphandle:
            kmphandle.write("reference,query,size,s,n,{}\n".format(','.join(header_names)))
            for query in query_set:
                for reference in reference_set:
                    sync_ibfjs = list()
                    for sampling_rate in args.sampling_rates:
                        query_name = kmp_dir.joinpath(kmp_wrapper.get_outname(query))
                        reference_name = kmp_dir.joinpath(kmp_wrapper.get_outname(reference))
                        sync_ibfj = kmp_wrapper.pairwise_jaccard(kmp_executable, query_name, reference_name)
                        sync_ibfjs.append(sync_ibfj)
                    kmphandle.write("{},{},{},{},{},{}\n".format(reference, query, args.sketch_size, s, n, ','.join(map(str, sync_ibfjs))))
    
    # [Post-processing] join results
    kmc_df = pd.read_csv(tmp_kmc_file, sep=",")
    mash_df = pd.read_csv(tmp_mash_file, sep=',')
    kmp_df = pd.read_csv(tmp_kmp_file, sep=',')
    joined = pd.merge(kmc_df, mash_df, on=["reference", "query"], how="inner")
    joined = pd.merge(joined, kmp_df, on=["reference", "query"], how="inner")
    joined = joined[["reference", "query", "size", "s", "n", "exj", "mhj"] + header_names] # drop unused columns
    joined.to_csv(args.output_csv, sep=',', header=True, index=False)

def extended_syncmers_experiment_main(args):
    kmp_executable = args.__bin_path.joinpath("kmp")
    cws_executable = args.__bin_path.joinpath("cws")
    assert kmp_executable.exists()
    assert 0 < args.k <= 32
    assert 0 <= args.z <= args.k
    assert args.extension > 0
    assert 3 <= args.repetitions <= 7
    assert (args.query_number != None and args.query_number >= 0) or (args.query_fraction != None and 0 <= args.query_fraction <= 1)
    assert args.sketch_size >= 0
    assert args.first.file() and args.second.is_file()
    assert args.tmp_dir.exists() and args.tmp_dir.is_dir()
    ek = args.extension + args.k
    assert ek <= 32

    wdir = args.tmp_dir
    extkmp_dir = wdir.joinpath("extended_syncmers_ibf")
    os.makedirs(extkmp_dir, exist_ok=True)
    n = kmp_wrapper.get_maximum_difference_from_byte_size(args.sketch_size, ek, args.repetitions, args.epsilon)

    random.seed(args.seed)
    tmp_i_sketch = wdir.joinpath(essentials.get_random_name("iblt.bin"))
    kmp_wrapper.sketch(kmp_executable, args.first, tmp_i_sketch, n, args.k, args.z, args.repetitions, args.epsilon, args.seed, wdir, args.max_ram)
    tmp_j_sketch = wdir.joinpath(essentials.get_random_name("iblt.bin"))
    kmp_wrapper.sketch(kmp_executable, args.second, tmp_j_sketch, n, args.k, args.z, args.repetitions, args.epsilon, args.seed, wdir, args.max_ram)
    #TODO
    

def minimizers_vs_syncmers_experiment_main(args):
    None

def main(args):
    if args.command == "jaccard": return jaccard_experiment_main(args)
    elif args.command == "extended": return extended_syncmers_experiment_main(args)
    elif args.command == "sampling": return minimizers_vs_syncmers_experiment_main(args)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_jaccard = subparsers.add_parser("jaccard", help="Compute multiple methods for estimating Jaccard")
    # parser_jaccard.add_argument("-a", "--algorithm", help="algorithms choices are: [exact (KMC), minHash (MASH), iblt]", type=str, choices=jaccard_choices, required=True)
    parser_jaccard.add_argument("-i", "--input-folder", help="folder containing fasta files to be sketched", type=pathlib.Path, required=True)
    parser_jaccard.add_argument("-o", "--output-csv", help="output csv file", type=pathlib.Path, required=True)
    parser_jaccard.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_jaccard.add_argument("-z", help="syncmer selection parameter when using syncmers", type=int, required=True)
    parser_jaccard.add_argument("-c", "--canonical", help="activate canonical k-mers", action="store_true")
    parser_jaccard.add_argument("-r", "--repetitions", help="number of repetitions in the IBLTs", type=int, default=4)
    parser_jaccard.add_argument("-e", "--epsilon", help="epsilon parameter of IBLTs", type=float, default=0)
    parser_jaccard.add_argument("-m", "--sketch-size", help="sketch size in Bytes for both minHashes and IBFs", type=int, required=True)
    query_option_parser = parser_jaccard.add_mutually_exclusive_group(required = True)
    query_option_parser.add_argument("-q", "--query-number", help="number of randomly selected sequences to be used as queries (disables -f option). 0 use all sequences (all-vs-all).", type=int)
    query_option_parser.add_argument("-f", "--query-fraction", help="fraction of input sequences to be considered as queries (disables -q option). 0 use all sequences (all-vs-all).", type=float)
    parser_jaccard.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)
    parser_jaccard.add_argument("-M", "--max-ram", help="maximum resident memory (GB) [4]", type=int, default=4)
    parser_jaccard.add_argument("-d", "--tmp-dir", help="folder where temporary files are stored", type=pathlib.Path, required=True)
    parser_jaccard.add_argument("--force-kmc", help="force recomputation of KMC databases", action="store_true")
    parser_jaccard.add_argument("--force-mash", help="force recomputation of mash output", action="store_true")
    parser_jaccard.add_argument("--force-kmp", help="force recomputation of kmp sketches and results", action="store_true")

    parser_extended = subparsers.add_parser("extended", help="Verify retrieval using extended syncmers")
    parser_extended.add_argument("-i", "--first", help="first fasta file", type=pathlib.Path, required=True)
    parser_extended.add_argument("-j", "--second", help="second fasta file", type=pathlib.Path, required=True)
    parser_extended.add_argument("-o", "--output-csv", help="output csv file to append to", type=pathlib.Path, required=True)
    parser_extended.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_extended.add_argument("-x", "--extension", help="number of bases added to the right of each k-mer", type=int, required=True)
    parser_extended.add_argument("-z", help="syncmer selection parameter when using syncmers", type=int, required=True)
    parser_extended.add_argument("-c", "--canonical", help="activate canonical k-mers", action="store_true")
    parser_extended.add_argument("-r", "--repetitions", help="number of repetitions in the IBLTs", type=int, default=4)
    parser_extended.add_argument("-e", "--epsilon", help="epsilon parameter of IBLTs", type=float, default=0)
    parser_extended.add_argument("-m", "--sketch-size", help="sketch size in Bytes for both minHashes and IBFs", type=int, required=True)
    parser_extended.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)
    parser_extended.add_argument("-M", "--max-ram", help="maximum resident memory (GB) [4]", type=int, default=4)
    parser_extended.add_argument("-d", "--tmp-dir", help="folder where temporary files are stored", type=pathlib.Path, required=True)

    return parser

if __name__ == "__main__":
    parser = parser_init()
    args = parser.parse_args(sys.argv)

    this = pathlib.Path(sys.argv[0])
    assert this.exists()
    this_path = this.parent.absolute()
    assert this_path.name == "scripts"
    repo_path = this_path.parent
    kmp_path = repo_path.joinpath("debug_build")
    assert kmp_path.exists()

    setattr(args, "__bin_path", kmp_path)
    setattr(args, "sampling_rates", [1])
    main(args)