from argparse import ArgumentError
import multiprocessing
from operator import sub
import os
import sys
import time
import math
import string
import random
import pathlib
import subprocess
from typing import TextIO
import pandas as pd
import fastx
import kmer
import mutate
import pbench

main_exec = "ibltseq"
toset_script = "2set.py"

kmc_cws_exec = "cws"

def get_random_name(extension: str, size=20, chars=string.ascii_uppercase + string.digits) -> str:
    suffix = ('.' + extension) if extension else ""
    return "".join(random.choice(chars) for _ in range(size)) + suffix

def txt2fasta(txt: str, fasta: str):
    counter = 0
    with open(txt, "r") as sh:
        with open(fasta, "w") as fh:
            for line in sh:
                fh.write(">{}\n{}\n".format(counter, line.strip()))#line has \n at the end already
                counter += 1

def isempty(file: str):
    try:
        with open(file, "r") as csvh:
            csvh.seek(0, os.SEEK_END)
            if not csvh.tell(): write_header = True
            else: write_header = False
    except FileNotFoundError as e:
        write_header = True
    return write_header

def get_path_leaf(path) -> tuple[str]:
    head, tail = os.path.split(path)
    if tail: return head, tail
    else: return head, os.path.basename(head)

def get_fastx_format(filename):
    if filename.endswith("fa"): return ("fasta", False)
    elif filename.endswith("fq"): return ("fastq", False)
    elif filename.endswith("fna"): return ("fasta", False)
    elif filename.endswith("fnq"): return ("fastq", False)
    elif filename.endswith("fasta"): return ("fasta", False)
    elif filename.endswith("fastq"): return ("fastq", False)
    elif filename.endswith("fa.gz"): return ("fasta", True)
    elif filename.endswith("fq.gz"): return ("fastq", True)
    elif filename.endswith("fna.gz"): return ("fasta", True)
    elif filename.endswith("fnq.gz"): return ("fasta", True)
    elif filename.endswith("fasta.gz"): return ("fasta", True)
    elif filename.endswith("fastq.gz"): return ("fasta", True)
    else: raise ValueError("Input file not in fastx format")

def get_fastx_in_folder(path: str) -> list[str]:
    fastx_files = list()
    for f in os.listdir(path):
        if os.path.isfile(os.path.join(path, f)):
            try:
                get_fastx_format(f)
                fastx_files.append(os.path.join(path, f))
            except ValueError:
                pass
    fastx_files.sort(key = lambda v : v.upper())
    return fastx_files

def get_kmc_paths(k: int, input_file: str, output_path: str):
    _, filename = get_path_leaf(input_file)
    filename = filename.split('.')[0]
    
    cwd = os.getcwd()
    kfolder = os.path.join(cwd, output_path)
    kfolder = os.path.join(kfolder, "k"+str(k))
    intermidiate_name = "dummy_" + filename
    intermidiate_file = os.path.join(kfolder, intermidiate_name)
    sorted_name = "S_" + filename
    sorted_file = os.path.join(kfolder, sorted_name)
    final_name = filename
    final_file = os.path.join(kfolder, final_name)
    return filename, kfolder, intermidiate_file, sorted_file, final_file

def kmc_count(k: int, input_file: str, output_path: str, kmc_dummy_folder: str, mmemory: int, unsorted: bool, canonical: bool) -> str:
    """Call kmc for k-mer counting"""
    _, kfolder, intermidiate_file, sorted_file, final_file = get_kmc_paths(k, input_file, output_path)
    fmt, _ = get_fastx_format(input_file)
    if fmt == "fasta": fmt = "-fm"
    elif fmt == "fastq": fmt = "-fq"
    else: raise ValueError("Unrecognized file type")
    try:
        os.makedirs(kfolder)
    except FileExistsError:
        pass
    this_path = os.path.dirname(os.path.abspath(__file__))
    # kmc = os.path.join(this_path, "kmc")
    # kmc_tools = os.path.join(this_path, "kmc_tools")
    kmc_command = ["kmc"]
    if not canonical: kmc_command.extend(["-b"])
    kmc_command.extend(["-k"+str(k), "-ci0", "-cs4294967295", "-cx4294967295", "-m"+str(mmemory), fmt, input_file, intermidiate_file, kmc_dummy_folder])
    out = subprocess.run(kmc_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if out.returncode != 0: exit(1)
    if (not unsorted): 
        out = subprocess.run(["kmc_tools", "transform", intermidiate_file, "sort", sorted_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if out.returncode != 0: exit(1)
    if os.path.exists(sorted_file+".kmc_pre"): os.remove(intermidiate_file+".kmc_pre")
    else: os.rename(intermidiate_file+".kmc_pre", sorted_file+".kmc_pre")
    if os.path.exists(sorted_file+".kmc_suf"): os.remove(intermidiate_file+".kmc_suf")
    else: os.rename(intermidiate_file+".kmc_suf", sorted_file+".kmc_suf")
    os.replace(sorted_file+".kmc_pre", final_file+".kmc_pre")
    os.replace(sorted_file+".kmc_suf", final_file+".kmc_suf")
    return final_file

def kmc_count_folder(infolder: str, outfolder: str, ks: list, kmc_dummy_folder: str, mmemory: int, unsorted: bool, canonical: bool) -> list[str]:
    """Build kmc databases for the input folder for all k values.
    
    The output folder will have a folder of the same name as the input with multiple folders named as k<value>.
    This function also creates the temporary folder for kmc hard-disk computations.
    """
    assert os.path.isdir(infolder)
    assert type(ks) == type(list())
    assert mmemory > 0
    fastxs = get_fastx_in_folder(infolder)
    _, basename = get_path_leaf(infolder)
    working_folder = os.path.join(outfolder, basename)
    try: os.makedirs(working_folder)
    except FileExistsError: pass
    try: os.makedirs(kmc_dummy_folder)
    except FileExistsError: pass
    kmced_files = list()
    for k in ks:
        for fpath in fastxs:
            kmced_files.append(kmc_count(k, fpath, working_folder, kmc_dummy_folder, mmemory, unsorted, canonical))

def kmc_jaccard(executable: str, kmc_file1: str, kmc_file2: str) -> tuple[float, float, float]:
    cws_command = [executable, "full", kmc_file1, kmc_file2]
    out = subprocess.run(cws_command, stderr=sys.stderr, stdout=subprocess.PIPE)
    if out.returncode != 0:
        sys.stderr.write("Error while performing kmc jaccard comparison\n")
        sys.exit(os.EX_SOFTWARE)
    return tuple(map(float, out.stdout.decode("utf-8").split(',')))

def mash_sketch(reference_set: list[str], output_path: str, k: int, s: int, seed: int, canonical: bool):
    # use mash sketch of multiple sequences + mash dist of multiple files
    sketch_command = ["mash", "sketch", "-o", output_path, "-k", str(k), "-s", str(s), "-S", str(seed)] #ignore strand
    if not canonical: sketch_command.extend(["-n"])
    sketch_command.extend(reference_set)
    out = subprocess.run(sketch_command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if out.returncode != 0:
        sys.stderr.write("Error, mash sketch command returned with exit code {}\n".format(out.returncode))
        sys.exit(os.EX_CANTCREAT)

def mash_dist(reference: str, queries: list, k: int, s: int, seed: int, tmp_file):
    dist_command = ["mash", "dist", "-s", str(s), "-S", str(seed), "-p", str(multiprocessing.cpu_count())]
    if not reference.endswith(".msh"): dist_command.extend(["-k", str(k)])
    dist_command.append(reference)
    dist_command.extend(queries)
    out = subprocess.run(dist_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if out.returncode != 0:
        sys.stderr.write("Error from mash dist command\n")
        sys.exit(os.EX_CANTCREAT)
    lines = out.stdout.decode("utf-8").split('\n')
    with open(tmp_file, "w") as mhandle:
        mhandle.write("reference,query,mhj\n")
        for line in lines:
            if line:
                rf, qr, _, _, frac = line.split('\t')
                n, d = frac.split('/')
                j = int(n)/int(d)
                mhandle.write("{},{},{}\n".format(rf,qr,j))

def configure4sequences(path: str, algorithm: str, k: int, z: int, store_hashes: bool, debug: bool):
    configure_path = os.path.dirname(path)
    configure_command = ["python3", os.path.join(configure_path, "configure.py"), "sequences", "-a", algorithm, "-k", str(k), "-z", str(z)] + [] if not store_hashes else ["--store_hashes"] + ["debug"] if debug else []
    out = subprocess.run(configure_command, stdout=subprocess.DEVNULL)
    if out.returncode != 0: 
        sys.stderr.write("Unable to configure ibltseq for sequence storage")
        exit(out.returncode)

def configure4variable_length_sequences(path: str, l: int, store_hashes: bool, debug: bool):
    configure_path = os.path.dirname(path)
    configure_command = ["python3", os.path.join(configure_path, "configure.py"), "vls", "-l", str(l)] + [] if not store_hashes else ["--store_hashes"] + ["debug"] if debug else []
    out = subprocess.run(configure_command, stdout=subprocess.DEVNULL)
    if out.returncode != 0: 
        sys.stderr.write("Unable to configure ibltseq for fragment storage")
        exit(out.returncode)

def configure4hashes(path: str, hash_width: int, minhash_dimension: int):
    configure_path = os.path.dirname(path)
    configure_command = ["python3", os.path.join(configure_path, "configure.py"), "hashes", "-l", str(hash_width), "-s", str(minhash_dimension)]
    out = subprocess.run(configure_command, stdout=subprocess.DEVNULL)
    if out.returncode != 0:
        sys.stderr.write("Unable to configure ibltseq for hash storage")
        exit(out.returncode)

def configure_length(algorithm: str, k: int, m: int):
    '''FIXME possible source of bugs beacause it is a copy of the same function in configure.py'''
    assert k > 0
    assert 0 <= m <= k
    if algorithm == "kmers": return k
    elif algorithm == "minimizers": return 2*k-m #space for k-mers grouped together
    elif algorithm == "syncmers": return k #allocate space to store syncmers
    elif algorithm == "segmentation": return 2*k #allocate enough space for the segmentation algorithm
    else: raise ValueError("This should never happen because the parser should catch it with options")

def get_ibfon(algorithm: str, sampling_rate: int, outpath: str, inputfile: str):
    bname = os.path.basename(inputfile)
    rname = bname.split('.')[0]
    return os.path.join(outpath, "{}_{}_r{}.ibf.bin".format(rname, algorithm, sampling_rate))

def jaccard_helper(aldiff_path: str, method: str, seq1: str, seq2: str, k: int, m: int, r: int, seed: int):
    executable = os.path.join(aldiff_path, main_exec)
    toset = os.path.join(aldiff_path, toset_script)
    dummy_fasta_header = ">0\n{}"
    if method == "minimizers": opt = 'w' #use windowed minimizers
    else: opt = 'm'
    if r == 1: command = "{} {} -k {} -{} {} -S {} | python3 {}".format(executable, method, k, opt, m, seed, toset)#sampling off
    else: command = "{} {} -k {} -{} {} -S {} | {} sample -r {} | python3 {}".format(executable, method, k, opt, m, seed, executable, r, toset)
    # sys.stderr.write("method: {}\n\n".format(method))
    # print(method)
    out = subprocess.run(command, text=True, input=dummy_fasta_header.format(seq1), stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    if out.returncode != 0: raise RuntimeError #Does not work because we have a pipe, so no way to check the output of executable
    s1 = set(out.stdout.split('\n'))
    if "" in s1: s1.remove("")
    # sys.stderr.write("table1:\n")
    # sys.stderr.write('\n'.join(sorted(list(s1))))
    # sys.stderr.write("\n")
    out = subprocess.run(command, text=True, input=dummy_fasta_header.format(seq2), stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    if out.returncode != 0: raise RuntimeError
    s2 = set(out.stdout.split('\n'))
    if "" in s2: s2.remove("")
    # sys.stderr.write("table2:\n")
    # sys.stderr.write('\n'.join(sorted(list(s2))))
    # sys.stderr.write("\n")
    intersection = len(s1 & s2)
    union = len(s1 | s2)
    # print("size (table1, table2) = ({}, {})".format(len(s1), len(s2)))
    # print("intersection/union = {}/{} = {}".format(intersection, union, intersection/union))
    # sys.stderr.write("\n---------------------------------------------\n")
    assert union
    return intersection / union

def ibf_sketch_set(executable: str, input_file: str, output_file: str, frag_len: int, n: int, r: int, epsilon: float, seed: int):
    assert input_file
    assert output_file
    assert 0 < frag_len <= 32
    assert n > 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    input_opt = ["-i", input_file]
    build = [executable, "build"]
    output_opt = ["-o", output_file]
    k_opt = ["-l", str(frag_len)]
    n_opt = ["-n", str(n)]
    r_opt = ["-r", str(r)]
    epsilon_opt = ["-e", str(epsilon)]
    seed_opt = ["-s", str(seed)]
    command = build + input_opt + output_opt + k_opt + n_opt + r_opt + epsilon_opt + seed_opt
    out = subprocess.run(command, stdout=sys.stdout, stderr=sys.stderr)
    if out.returncode != 0: sys.stderr.write("ibltseq exited with error {}\n".format(out.returncode))
    
def ibf_sketch_each_sequence(executable: str, algo: str, input_fastx: str, output_folder: str, k: int, m: str, n: int, r: int, epsilon: float, seed: int):
    #TODO it needs a little bit of refactoring
    assert input_fastx
    assert output_folder and os.path.isdir(output_folder)
    assert 0 < k <= 32
    assert 0 <= m < k
    assert n > 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    aldiff_path = os.path.dirname(executable)
    build = [executable, "build"]
    k_opt = ["-k", str(k)]
    n_opt = ["-n", str(n)]
    r_opt = ["-r", str(r)]
    epsilon_opt = ["-e", str(epsilon)]
    seed_opt = ["-s", str(seed)]
    partial_build_command = build + k_opt + n_opt + r_opt + epsilon_opt + seed_opt
    input_basename = pathlib.Path(input_fastx).stem
    dummy_fasta_name = get_random_name("fna")
    dummy_input_name = get_random_name("txt")
    with open(dummy_fasta_name, "w") as dh:
        with open(input_fastx, "r") as fh:
            for name, seq, _ in fastx.read(fh):
                dh.seek(0)
                dh.write(">{}\n{}".format(name, seq))
                dh.truncate()
                dh.flush()
                fragment_command = "{} {} -i {} -k {} -m {} | python3 {} > {}".format(executable, algo, dummy_fasta_name, k, m, os.path.join(aldiff_path, toset_script), dummy_input_name)
                out = subprocess.run(fragment_command, capture_output=False, shell=True)
                if out.returncode != 0: 
                    sys.stderr.write("Unable to compute {} table for record {} (failed with exit code {})\n".format(algo, name, out.returncode))
                else:
                    output_name = "{}_{}.{}".format(input_basename, name, "ibf.bin")
                    output_opt = ["-o", os.path.join(output_folder, output_name)]
                    input_opt = ["-i", dummy_input_name]
                    out = subprocess.run(partial_build_command + input_opt + output_opt, stdout=sys.stdout, stderr=sys.stderr)
                    if out.returncode != 0: sys.stderr.write("ibltseq construction for record {} failed with exit code {}\n", name, out.returncode)

# def ibf_sketch_fastx(executable: str, toset: str, algorithm: str, input_file: str, output_file: str, k: int, z: int, n: int, r: int, epsilon: float, seed: int, sampling_rate: int, unused: int, canonical: bool):
#     """WORKING VERSION"""
#     assert not (extra_len and canonical) #having canonical form of groups is a waste of time because k-mers must be extracted later anyway
#     fragmentation_command = [executable]
#     if algorithm == "syncmers": fragmentation_command.extend(["syncmers", "-m", str(z)])
#     elif algorithm == "minimizers": fragmentation_command.extend(["minimizers", "-m", str(z)])
#     elif algorithm == "kmers": fragmentation_command.extend(["kmers"])
#     else: pass #This should never happen because configure_length raise error if algorithm unrecognized
#     # if sampling_rate != 1: #skip sample command if no sampling
#     #     sampling_command = [executable, "sample", "-r", str(sampling_rate)]
#     fragmentation_command.extend(["-i", input_file, "-k", str(k)])
#     toset_command = ["python3", toset]
#     build_command = [executable, "build", "-o", output_file, "-n", str(n), "-r", str(r), "-e", str(epsilon), "-s", str(seed)]
#     # if (canonical): toset_command.extend(["-c"])
#     if canonical:
#         very_bad_workaround = "{} kmers -i {} -k {} | python3 {} -c | python3 kmer.py syncmers -z {} | {} build -o {} -n {} -r {} -e {} -s {}".format(executable, input_file, k, toset, z, executable, output_file, n, r, epsilon, seed)
#         build_out = subprocess.run(very_bad_workaround, shell=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
#     else:
#         fragmentation_out = subprocess.Popen(fragmentation_command, stdout=subprocess.PIPE, stderr=sys.stderr)
#         if sampling_rate != 1: 
#             sampling_command = [executable, "sample", "-r", str(sampling_rate)]
#             sampling_out = subprocess.Popen(sampling_command, stdin=fragmentation_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
#         else: sampling_out = fragmentation_out
#         toset_out = subprocess.Popen(toset_command, stdin=sampling_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
#         build_out = subprocess.Popen(build_command, stdin=toset_out.stdout, stdout=sys.stdout, stderr=sys.stderr)
#         # command = fragmentation_command + ['|'] + sampling_command + (['|'] if sampling_rate != 1 else []) + toset_command + ['|'] + build_command
#         # out = subprocess.run(command, stdout=sys.stdout, stderr=sys.stderr, shell=True)
#         build_out.communicate()
#         toset_out.communicate()
#         fragmentation_out.communicate()
#     if build_out.returncode != 0:
#         sys.stderr.write("building from fastx file exited with error\n")

def fasta_to_token_fasta(executable: str, toset: str, algorithm: str, input_file: str, output_file: str, k: int, z: int, sampling_rate: int, canonical: bool):
    fragmentation_command = [executable]
    if algorithm == "syncmers" and not canonical: fragmentation_command.extend(["syncmers", "-m", str(z)])
    elif algorithm == "minimizers": fragmentation_command.extend(["minimizers", "-m", str(z)])
    elif canonical or algorithm == "kmers": fragmentation_command.extend(["kmers"])
    else: sys.stderr.write("This should never happen\n") 
    fragmentation_command.extend(["-i", input_file, "-k", str(k)])
    toset_command = ["python3", toset]
    tofasta_command = ["python3", "kmer.py", "2fasta", "-o", output_file]
    fragmentation_out = subprocess.Popen(fragmentation_command, stdout=subprocess.PIPE, stderr=sys.stderr)
    if canonical: #toset_command.append("-c")
        canonical_out = subprocess.Popen(["python3", toset, "-c"], stdin=fragmentation_out.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        syncmer_sampling_out = subprocess.Popen(["python3", "kmer.py", "syncmers", "-z", str(z)], stdin=canonical_out.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    else:
        syncmer_sampling_out = fragmentation_out
    if sampling_rate != 1: 
        sampling_command = [executable, "sample", "-r", str(sampling_rate)]
        sampling_out = subprocess.Popen(sampling_command, stdin=syncmer_sampling_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    else: sampling_out = syncmer_sampling_out
    toset_out = subprocess.Popen(toset_command, stdin=sampling_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    tofasta_out = subprocess.Popen(tofasta_command, stdin=toset_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    tofasta_out.communicate()
    if tofasta_out.returncode != 0:
        sys.stderr.write("conversion of fasta to sampled fasta exited with error\n")

def ibf_sketch_fastx(executable: str, toset: str, algorithm: str, input_file: str, output_file: str, k: int, z: int, n: int, r: int, epsilon: float, seed: int, sampling_rate: int, extra_len: int, canonical: bool):
    """FINAL VERSION"""
    alternative_syncmer_sampling = None
    fragmentation_command = [executable]
    if algorithm == "syncmers" and not canonical: fragmentation_command.extend(["syncmers", "-m", str(z), "-g", str(extra_len)])
    elif algorithm == "minimizers": fragmentation_command.extend(["minimizers", "-m", str(z)])
    elif algorithm == "kmers" or (canonical and not extra_len) : fragmentation_command.extend(["kmers"])
    elif canonical and extra_len: raise ValueError("canonical extended syncmers are not possible here")
    else: sys.stderr.write("This should never happen\n") #This should never happen because configure_length raise error if algorithm unrecognized
    fragmentation_command.extend(["-i", input_file, "-k", str(k)])
    toset_command = ["python3", toset]
    build_command = [executable, "build", "-o", output_file, "-n", str(n), "-r", str(r), "-e", str(epsilon), "-s", str(seed)]
    fragmentation_out = subprocess.Popen(fragmentation_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if canonical: #toset_command.append("-c")
        # Computing syncmers and then taking their cononical form doesn't work because we miss all syncmers in canonical form given the sampling
        # The solution is to compute the whole set of canonical k-mers and then sample it using syncmers
        # Here we use python scripts because the syncmer subcommand of ibltseq would need to do something similar so better use what we already have.
        canonical_out = subprocess.Popen(["python3", toset, "-c"], stdin=fragmentation_out.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        syncmer_sampling_out = subprocess.Popen(["python3", "kmer.py", "syncmers", "-z", str(z)], stdin=canonical_out.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        #very_bad_workaround = "{} kmers -i {} -k {} | python3 {} -c | python3 kmer.py syncmers -z {} | {} build -o {} -n {} -r {} -e {} -s {}".format(executable, input_file, k, toset, z, executable, output_file, n, r, epsilon, seed)
        #build_out = subprocess.run(very_bad_workaround, shell=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    else:
        syncmer_sampling_out = fragmentation_out
    if sampling_rate != 1: 
        sampling_command = [executable, "sample", "-r", str(sampling_rate)]
        sampling_out = subprocess.Popen(sampling_command, stdin=syncmer_sampling_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    else: sampling_out = syncmer_sampling_out
    toset_out = subprocess.Popen(toset_command, stdin=sampling_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    build_out = subprocess.Popen(build_command, stdin=toset_out.stdout, stdout=sys.stdout, stderr=sys.stderr)
    build_out.communicate()
    toset_out.communicate()#just to be sure
    if build_out.returncode != 0:
        sys.stderr.write("building from fastx file exited with error\n")

def compute_pairwise_jaccard_from_ibfs(executable: str, sk1: str, sk2: str) -> float:
    out = subprocess.run([executable, "jaccard", "-i", sk1, "-j", sk2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out.returncode != 0: return float("NaN")
    #   sys.stderr.write(".")
    #   sys.stderr.write("Unable to compute IBF Jaccard between:\n {}\n{}\n".format(sk1, sk2))
    # exit(1)
    return float(out.stdout.decode("utf-8").split(',')[0])

def compute_jaccard_from_ibfs(executable: str, sketch_list: list, rh: TextIO):
    #NOTE this function does not use compute_pairwise_jaccard_from_ibfs for optimization purposes
    assert sketch_list
    assert rh
    n = len(sketch_list)
    partial_command = [executable, "jaccard"]
    for i in range(n):
        first_opt = ["-i", sketch_list[i]]
        name1 = os.path.basename(sketch_list[i]).rsplit('.', 2)[0]
        for j in range(i + 1, n):
            second_opt = ["-j", sketch_list[j]]
            out = subprocess.run(partial_command + first_opt + second_opt, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            # if out.returncode != 0: 
            #     sys.stderr.write("ibltseq jaccard computation between {} and {} failed with exit code {}\n".format(sketch_list[i], sketch_list[j], out.returncode))
            output = out.stdout.decode("utf-8") #This should be a csv line terminating with a \n
            name2 = os.path.basename(sketch_list[j]).rsplit('.', 2)[0]
            rh.write("{},{},{}".format(name1, name2, output))

def sketch_diff(executable: str, ibfi: str, ibfj: str, diff: str):
    diff_command = [executable, "diff", "-i", ibfi, "-j", ibfj, "-o", diff]
    out = subprocess.run(diff_command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if out.returncode != 0: 
        sys.stderr.write("Error from diff command\n")
        sys.exit(os.EX_CANTCREAT)

def sketch_list(executable: str, diff: str) -> bytes:
    list_command = [executable, "list", "-i", diff]
    out = subprocess.run(list_command, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
    if out.returncode != 0: 
        sys.stderr.write("Error from list command\n")
        sys.exit(os.EX_CANTCREAT)
    return out.stderr

def check_peelability(message: bytes) -> bool:
    return not (message.decode("utf-8") == "Warning: unpeelable sketch\n")

def compare_ibf_methods(aldiff_path: str, original: str, mutated: str, tmp_names: tuple[str, str, str, str, str, str, str], k: int, m: int, starting_n: int, max_trials: int) -> tuple[int, int]:
    '''Mutate sequence and test IBFs for the difference'''
    # NOTE ibf_sketch_fastx is explicitly unrolled here because we receive a sequence as input, not a fastx file
    executable = os.path.join(aldiff_path, main_exec)
    toset = os.path.join(aldiff_path, toset_script)
    dummy_mm_ori_table_name = tmp_names[0]
    dummy_sm_ori_table_name = tmp_names[1]
    dummy_mm_mut_table_name = tmp_names[2]
    dummy_sm_mut_table_name = tmp_names[3]
    dummy_ori_ibf = tmp_names[4]
    dummy_mut_ibf = tmp_names[5]
    dummy_tmp_ibf = tmp_names[6]
    minimizers_command = "{} minimizers -k {} -m {} | python3 {} > ".format(executable, k, m, toset)
    syncmers_command =   "{} syncmers   -k {} -m {} | python3 {} > ".format(executable, k, m, toset)
    out = subprocess.run(minimizers_command + dummy_mm_ori_table_name, capture_output=False, shell=True, input=">0\n"+original, text=True)
    if out.returncode != 0: 
        sys.stderr.write("Error while generating minimizers table for the original sequence\n")
        sys.exit(os.EX_CANTCREAT)
    out = subprocess.run(syncmers_command + dummy_sm_ori_table_name, capture_output=False, shell=True, input=">0\n"+original, text=True)
    if out.returncode != 0: 
        sys.stderr.write("Error while generating syncmers table for the original sequence\n")
        sys.exit(os.EX_CANTCREAT)
    out = subprocess.run(minimizers_command + dummy_mm_mut_table_name, capture_output=False, shell=True, input=">0\n"+mutated, text=True)
    if out.returncode != 0: 
        sys.stderr.write("Error while generating minimizers table for the mutated sequence\n")
        sys.exit(os.EX_CANTCREAT)
    out = subprocess.run(syncmers_command + dummy_sm_mut_table_name, capture_output=False, shell=True, input=">0\n"+mutated, text=True)
    if out.returncode != 0:
        sys.stderr.write("Error while generating syncmers table for the mutated sequence\n")
        sys.exit(os.EX_CANTCREAT)
    peelable = False
    n = starting_n
    mm_len = configure_length("minimizers", k, m)
    while not peelable and n < max_trials:
        ibf_sketch_set(executable, dummy_mm_ori_table_name, dummy_ori_ibf, mm_len, n, 3, 0, 42)
        ibf_sketch_set(executable, dummy_mm_mut_table_name, dummy_mut_ibf, mm_len, n, 3, 0, 42)
        sketch_diff(executable, dummy_ori_ibf, dummy_mut_ibf, dummy_tmp_ibf)
        message = sketch_list(executable, dummy_tmp_ibf)
        peelable = check_peelability(message)
        n += 1
    mm_n = n-1
    peelable = False
    n = starting_n
    sync_len = configure_length("syncmers", k, m)
    while not peelable and n < max_trials:
        ibf_sketch_set(executable, dummy_sm_ori_table_name, dummy_ori_ibf, sync_len, n, 3, 0, 42)
        ibf_sketch_set(executable, dummy_sm_mut_table_name, dummy_mut_ibf, sync_len, n, 3, 0, 42)
        sketch_diff(executable, dummy_ori_ibf, dummy_mut_ibf, dummy_tmp_ibf)
        message = sketch_list(executable, dummy_tmp_ibf)
        peelable = check_peelability(message)
        n += 1
    sm_n = n-1
    return mm_n, sm_n

def ibf_build_minhash_collection(executable: str, hash_width: int, minhash_size: int, n: int, alice: list[str], bob: list[str], output_file: str):
    command = [executable, "collection", "-w", str(hash_width), "-z", str(minhash_size), "-n", str(n), "-o", output_file, "-a"] + alice + ["-b"] + bob
    out = subprocess.run(command, stderr=sys.stderr, stdout=subprocess.PIPE)
    if out.returncode != 0: 
        sys.stderr.write("Error from collection command\n")
        sys.exit(os.EX_CANTCREAT)
    return int(out.stdout.decode("utf-8"))

def benchmark(command: list, special_code: int = 0) -> tuple[int, int, int, bytes]:
    if special_code < 0 or special_code > 255: sys.stderr.write("[Benchmark warning]: The list of return codes to be ignored contains invalid codes \n")
    bench = pbench.ProcessBenchmark(command)
    try:
        bench.execute()
        while bench.poll():
            time.sleep(0.000001)
    except OSError:
        sys.stderr.write("Unable to launch command: {}\n".format(' '.join(command)))
        sys.exit(1)
    else:
        _, err = bench.communicate()
    finally:
        bench.close()
    if bench.p.returncode == 0:
        return bench.t1 - bench.t0, bench.max_rss_memory, bench.max_vms_memory, True
    elif bench.p.returncode != 0 and bench.p.returncode == special_code:
        return bench.t1 - bench.t0, bench.max_rss_memory, bench.max_vms_memory, False
    else:
        sys.stderr.write("Error during benchmark (errorcode = {})\n".format(bench.p.returncode))
        sys.stderr.write("{}\n{}\n".format(command, err))
        sys.exit(os.EX_SOFTWARE)

def measure_ibf_method(aldiff_path: str, original: str, mutated: str, tmp_names: tuple[str, str, str, str, str, str, str, str], k: int, m: int, nestimated: int, epsilon: float, algorithm: str) -> tuple[int, int, int, int, int, int, int]:
    '''IBF statistics'''
    executable = os.path.join(aldiff_path, main_exec)
    toset = os.path.join(aldiff_path, toset_script)
    dummy_ori_table_name = tmp_names[0]
    dummy_mut_table_name = tmp_names[1]
    dummy_ori_ibf = tmp_names[2]
    dummy_mut_ibf = tmp_names[3]
    dummy_tmp_ibf = tmp_names[4]
    command = "{} {} -k {} -m {} | python3 {} > ".format(executable, algorithm, k, m, toset)
    frag_len = configure_length(algorithm, k, m)
    out = subprocess.run(command + dummy_ori_table_name, capture_output=False, shell=True, input=">0\n"+original, text=True)
    if out.returncode != 0: 
        sys.stderr.write("Error while generating the table for the original sequence\n")
        sys.exit(os.EX_CANTCREAT)
    out = subprocess.run(command + dummy_mut_table_name, capture_output=False, shell=True, input=">0\n"+mutated, text=True)
    if out.returncode != 0: 
        sys.stderr.write("Error while generating the table for the mutated sequence\n")
        sys.exit(os.EX_CANTCREAT)
    out = subprocess.run(["python3", "kmer.py", "tabdiff", "-i", dummy_ori_table_name, "-j", dummy_mut_table_name, "-s"], stdout=subprocess.PIPE, stderr=sys.stderr)
    if out.returncode != 0:
        sys.stderr.write("Error while computing true table differences\n")
        sys.exit(os.EX_CANTCREAT)
    ntrue = len(out.stdout.decode("utf-8").split('\n')) - 1
    build_command = [executable, "build", "-i", dummy_ori_table_name, "-o", dummy_ori_ibf, "-l", str(frag_len), "-n", str(nestimated), "-r", str(3), "-e", str(epsilon), "-s", str(42)]
    ori_build_time_ns, ori_build_rss, ori_build_vms, _ = benchmark(build_command)
    build_command[3] = dummy_mut_table_name
    build_command[5] = dummy_mut_ibf
    mut_build_time_ns, mut_build_rss, mut_build_vms, _ = benchmark(build_command)
    diff_command = [executable, "diff", "-i", dummy_ori_ibf, "-j", dummy_mut_ibf, "-o", dummy_tmp_ibf]
    diff_time_ns, diff_rss, diff_vms, _ = benchmark(diff_command)
    list_command = [executable, "list", "-i", dummy_tmp_ibf]
    list_time_ns, list_rss, list_vms, peelable = benchmark(list_command, 1) #1 because aldiff main hides all internal errors behind EXIT_FAILURE
    ibf_size = os.stat(dummy_tmp_ibf).st_size#in bytes
    # peelable = check_peelability(message)
    return ori_build_time_ns, ori_build_rss, ori_build_vms, mut_build_time_ns, mut_build_rss, mut_build_vms, diff_time_ns, diff_rss, diff_vms, list_time_ns, list_rss, list_vms, ibf_size, peelable, nestimated, ntrue

def measure_kmc_method(l: int, tmp_names: tuple[str, str, str, str, str, str, str, str]):
    # executable = os.path.join(aldiff_path, kmc)
    dummy_ori_table_name = tmp_names[0]
    dummy_mut_table_name = tmp_names[1]
    dummy_fasta_name = tmp_names[5]
    dummy_kmc_output = tmp_names[6]
    kmc_tmp_folder = tmp_names[7]
    kmc_command = ["kmc", "-k{}".format(l), "-fm", "-b", "-r", "-m4", "-ci0", "-cs0"]
    txt2fasta(dummy_ori_table_name, dummy_fasta_name)
    kmc_ori_time_ns, kmc_ori_rss, kmc_ori_vms, _ = benchmark(kmc_command + [dummy_fasta_name, dummy_kmc_output, kmc_tmp_folder])
    kmc_ori_size = os.stat(dummy_kmc_output + ".kmc_pre").st_size + os.stat(dummy_kmc_output + ".kmc_suf").st_size
    txt2fasta(dummy_mut_table_name, dummy_fasta_name)
    kmc_mut_time_ns, kmc_mut_rss, kmc_mut_vms, _ = benchmark(kmc_command + [dummy_fasta_name, dummy_kmc_output, kmc_tmp_folder])
    kmc_mut_size = os.stat(dummy_kmc_output + ".kmc_pre").st_size + os.stat(dummy_kmc_output + ".kmc_suf").st_size
    return kmc_ori_time_ns, kmc_ori_rss, kmc_ori_vms, kmc_ori_size, kmc_mut_time_ns, kmc_mut_rss, kmc_mut_vms, kmc_mut_size

def experiment1(scripts_path: str, slen: int, k: int, m: int, mut_rate: float, indel_fraction: float, ext_prob: float, max_trials: int, dummy_names: list) -> tuple[int, int, int, int, int, int]:
    aldiff_path = os.path.dirname(scripts_path)
    seq = ''.join(random.choice('ACTG') for _ in range(slen))
    mseq, nsub, nins, ndel, tidl = mutate.mutate(seq, mut_rate, indel_fraction, ext_prob)
    mm_maxn, sm_maxn = compare_ibf_methods(aldiff_path, seq, mseq, dummy_names, k, m, nsub, max_trials)
    return nsub, nins, ndel, tidl, mm_maxn, sm_maxn

def experiment2(scripts_path: str, slen: int, k: int, m: int, epsilon: float, mut_rate: float, indel_fraction: float, ext_prob: float, algorithm: str, dummy_names: list) -> tuple[int, int, int, int, int, int, int]:
    aldiff_path = os.path.dirname(scripts_path)
    seq = ''.join(random.choice('ACTG') for _ in range(slen))
    mseq, nsub, nins, ndel, tidl = mutate.mutate(seq, mut_rate, indel_fraction, ext_prob)
    # multiplier = 1.9
    ccs = (k-m+1)/2
    coverage_rate = k/ccs
    estimated_n = max(int(coverage_rate * nsub * 2), 10) #each mutation is covered by coverage_rate sampled k-mers and 2 because of symmetry
    # print("{}*{}*{}={}".format(coverage_rate, nsub, multiplier, estimated_n))
    ibf_blob = measure_ibf_method(aldiff_path, seq, mseq, dummy_names, k, m, estimated_n, epsilon, algorithm)
    kmc_blob = measure_kmc_method(configure_length(algorithm, k, m), dummy_names)
    return [nsub, nins, ndel, tidl] + list(ibf_blob) + list(kmc_blob)

def experiment5(scripts_path: str, seed: int, k: int, m: int):#, seq1: str, seq2: str):
    aldiff_path = os.path.dirname(scripts_path)
    seq1 = ''.join(random.choice('ACTG') for _ in range(args.l))
    seq2, _, _, _, _ = mutate.mutate(seq1, args.mutationp, args.indelf, args.extensionp)
    # compute minimizer window length to obtain the same density as syncmers
    #ccs = (k-m+1)/2
    #cmm = (w+1)/2
    #css = cmm -> (k-m+1)/2 = (w+1)/2 -> w = k-m
    w = k-m
    sampling_rate = int((k-m+1)/2)
    truejacc = jaccard_helper(aldiff_path, "kmers", seq1, seq2, k, m, 1, seed)
    sampledjacc = jaccard_helper(aldiff_path, "kmers", seq1, seq2, k, m, sampling_rate, seed)
    mmjacc = jaccard_helper(aldiff_path, "minimizers", seq1, seq2, k, w, 1, seed)
    syncjacc = jaccard_helper(aldiff_path, "syncmers", seq1, seq2, k, m, 1, seed)
    return truejacc, sampledjacc, mmjacc, syncjacc

#--------------------------------------------------------------------------------------------------------------------------

def sketch_main(args):
    executable = os.path.join(os.path.dirname(args.__path), main_exec)
    ibf_sketch_set(executable, args.i, args.o, args.l, args.n, args.repetitions, args.epsilon, args.seed)

def multi_main(args):
    executable = os.path.join(os.path.dirname(args.__path), main_exec)
    ibf_sketch_each_sequence(executable, args.algorithm, args.i, args.o, args.k, args.m, args.n, args.repetitions, args.epsilon, args.seed)

def jaccard_main(args):
    executable = os.path.join(os.path.dirname(args.__path), main_exec)
    if len(args.i) == 1:
        mypath = args.i[0]
        args.i = [os.path.join(mypath, f) for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    else:
        args.i = [f for f in args.i if os.path.isfile(f)]
    with open(args.o, "w") as rh:
        compute_jaccard_from_ibfs(executable, args.i, rh)

def rndset_main(args):
    assert args.l >= 0
    assert args.n >= 0
    assert 0 <= args.mutationp <= 1
    assert 0 <= args.indelf <= 1
    assert 0 <= args.extensionp <= 1
    folder = os.path.isdir(args.o)
    random.seed(args.seed)
    seq = ''.join(random.choice('ACTG') for _ in range(args.l))
    if folder: fname = os.path.join(args.o, "original.fasta")
    else: fname = "{}_original.fasta".format(args.o, idx)
    with open(fname, "w") as fh:
        fastx.fasta_write(fh, "origin", seq, True)
    for idx in range(args.n):
        mseq, _, _, _, _ = mutate.mutate(seq, args.mutationp, args.indelf, args.extensionp)
        if folder: fname = os.path.join(args.o, "{}.fasta".format(idx))
        else: fname = "{}_{}.fasta".format(args.o, idx)
        with open(fname, "w") as fh:
            fastx.fasta_write(fh, str(idx), mseq, True)

def experiment1_main(args):
    '''Comparison between minimizers and syncmers on random sequences

    NOTE Each IBF cell has enough space for minimizer grouping, which works for syncmers too (2*k-m is always >= than k).
    For space estimation, see experiment2.
    Given the same initial seed, ex1 and ex2 must guarantee that the series of random genomes is the same.
    '''
    assert args.l >= 0
    assert 0 < args.k <= 32
    assert 0 <= args.z < args.k
    assert 0 <= args.mutationp <= 1
    assert 0 <= args.indelf <= 1
    assert 0 <= args.extensionp <= 1
    if not args.o: oh = sys.stdout
    else: oh = open(args.o, "w")
    configure4sequences(args.__exepath, "minimizers", args.k, args.z, False, False)
    random.seed(args.seed)
    mm_otn = get_random_name("txt")
    sm_otn = get_random_name("txt")
    mm_mtn = get_random_name("txt")
    sm_mtn = get_random_name("txt")
    ibf_on = get_random_name("ibf.bin")
    ibf_mn = get_random_name("ibf.bin")
    ibf_diff = get_random_name("ibf.bin")
    dummy_names = [mm_otn, sm_otn, mm_mtn, sm_mtn, ibf_on, ibf_mn, ibf_diff]
    # print(dummy_names)
    oh.write("substitutions,insertions,deletions,indellen,nmm,nsm\n") #CSV header
    for _ in range(args.n):
        pack = experiment1(args.__exepath, args.l, args.k, args.z, args.mutationp, args.indelf, args.extensionp, args.mtrials, dummy_names)
        oh.write("{}\n".format(','.join(map(str, pack))))
    oh.close()

def experiment2_main(args):
    '''Measure (time, RAM) X (construction, listing) and disk space on random sequences using IBFs and KMC.'''
    assert args.l >= 0
    assert 0 < args.k <= 32
    assert 0 <= args.z < args.k
    assert args.e >= 0
    assert 0 <= args.mutationp <= 1
    assert 0 <= args.indelf <= 1
    assert 0 <= args.extensionp <= 1
    configure4sequences(args.__exepath, args.algorithm, args.k, args.z, False, False)
    random.seed(args.seed)
    otn = os.path.join(args.wfolder, get_random_name("txt"))
    mtn = os.path.join(args.wfolder, get_random_name("txt"))
    ibf_on = os.path.join(args.wfolder, get_random_name("ibf.bin"))
    ibf_mn = os.path.join(args.wfolder, get_random_name("ibf.bin"))
    ibf_diff = os.path.join(args.wfolder, get_random_name("ibf.bin"))
    dummy_fasta = os.path.join(args.wfolder, get_random_name("fna"))
    dummy_kmc_output = os.path.join(args.wfolder, get_random_name(""))
    dummy_names = [otn, mtn, ibf_on, ibf_mn, ibf_diff, dummy_fasta, dummy_kmc_output, args.wfolder]
    if not args.o: oh = sys.stdout
    else: oh = open(args.o, "w")
    oh.write("length,k,z,epsilon,mutp,indelf,extp,substitutions,insertions,deletions,indellen,obtime,obrss,obvms,mbtime,mbrss,mbvms,dtime,drss,dvms,ltime,lrss,lvms,ibfsize,peelable,estimated n,true n,kmcotime,kmcorss,kmcovms,kmcosize,kmcmtime,kmcmrss,kmcmvms,kmcmsize\n")
    for _ in range(args.n):
        pack = experiment2(args.__exepath, args.l, args.k, args.z, args.e, args.mutationp, args.indelf, args.extensionp, args.algorithm, dummy_names)
        pack = [args.l, args.k, args.z, args.e, args.mutationp, args.indelf, args.extensionp] + pack
        oh.write("{}\n".format(','.join(map(str, pack))))
    oh.close()

def experiment3_main(args):
    '''Compute Jaccard using multiple techniques: full set comparisons (use mash with very big sketch size), minHash, naive sampling + IBFs, syncmers + IBFs.

    Use spearman correlation to compare: true values, values obtained by minHash sampling, values obtained using IBFs + sampling.
    How well do syncmers + IBF approximate the true mutation rate (under a simple poisson process)?
    '''
    assert 0 < args.k <= 32
    assert 0 <= args.z <= args.k
    assert (args.qnumber and args.qnumber >= 0) or (args.qfraction and 0 <= args.qfraction <= 1)
    assert args.ssize >= 0
    assert os.path.isdir(args.ifolder)
    aldiff_path = os.path.dirname(args.__exepath)
    ibltseq_exec = os.path.join(aldiff_path, main_exec)
    toset_exec = os.path.join(aldiff_path, toset_script)
    cws_exec = os.path.join(aldiff_path, kmc_cws_exec)

    random.seed(args.seed)
    tmp_mash_file = os.path.join(args.wfolder, get_random_name("csv"))
    tmp_other_file = os.path.join(args.wfolder, get_random_name("csv"))
    kmc_dir = os.path.join(args.wfolder, "kmc")
    _, dataset_basename = get_path_leaf(args.ifolder)
    kmc_dataset_dir = os.path.join(kmc_dir, dataset_basename) # used later to retrieve kmc databases
    kmc_tmp_dir = os.path.join(args.wfolder, "kmc_tmp")
    mash_dir = os.path.join(args.wfolder, "mash")
    mash_ref = os.path.join(mash_dir, get_random_name(""))
    tmp_syncmers_mash_file = os.path.join(args.wfolder, get_random_name("csv"))
    syncmers_mash_ref = os.path.join(mash_dir, get_random_name(""))
    syncmers_fasta_dir = os.path.join(args.wfolder, "syncmers_fastas")
    sync_ibf_dir = os.path.join(args.wfolder, "syncmer_ibf")
    extended_ibf_dir = os.path.join(args.wfolder, "extended_syncmers_ibf")
    os.makedirs(kmc_dir, exist_ok=True)
    os.makedirs(kmc_tmp_dir, exist_ok=True)
    os.makedirs(mash_dir, exist_ok=True)
    os.makedirs(sync_ibf_dir, exist_ok=True)
    os.makedirs(extended_ibf_dir, exist_ok=True)
    os.makedirs(syncmers_fasta_dir, exist_ok=True)

    mh_s_param = math.ceil(args.ssize / (4 if 4**args.k < 2**32 else 8)) # number of hashes in a minHash sketch (hashes of 8 bytes)
    payload_size = int((args.k * 2 + 8 - 1) / 8)
    counter_size = 8
    bucket_size = payload_size + counter_size
    repetitions = 3
    epsilon = 0
    peelability_oversize = 1.222 + epsilon
    ibf_n_param = math.ceil((args.ssize - bucket_size * repetitions) / (bucket_size * peelability_oversize)) # maximum number of differences trackable by the given sketch size
    # sampling_rate = math.ceil(args.L0 / ibf_n_param)
    # z = args.k + 1 - 2 * sampling_rate
    if mh_s_param == 0: raise RuntimeError("minHash sketches of size 0 given the current parameters")
    if ibf_n_param == 0: raise RuntimeError("IBF sketches of size 0 given the current parameters")
    # if z > args.k or z < 0: raise RuntimeError("The selected sketch space is too small to give meaningful syncmer sampling: z={}".format(z))
    # partition file set into query (Q) and reference (R)
    data_set = get_fastx_in_folder(args.ifolder)
    # if len(data_set) <= 1: sys.stderr.write("Too few input sequences to have meaningful results\n") # Write something to the csv?
    if args.qnumber and not args.qfraction: query_size = args.qnumber
    elif not args.qnumber and args.qfraction: query_size = max(1, int(len(data_set) * args.qfraction))
    else: raise RuntimeError("Something went wrong with the mutually-exclusive options 'qnumber' and 'qfraction'")
    query_set = sorted(random.sample(data_set, query_size))
    if (query_size == len(data_set)): reference_set = data_set
    else: reference_set = sorted(list(set(data_set) - set(query_set)))
    # Count files using KMC
    _ = kmc_count_folder(args.ifolder, kmc_dir, [args.k], kmc_tmp_dir, 8, False, args.canonical)
    # Build MASH sketches on input sequences
    mash_sketch(reference_set, mash_ref, args.k, mh_s_param, args.seed, args.canonical)
    mash_dist(mash_ref+".msh", query_set, args.k, mh_s_param, args.seed, tmp_mash_file)
    # Build MASH sketches on syncmers sets
    get_sampled_fasta = lambda f: os.path.join(syncmers_fasta_dir, os.path.basename(f))
    for f in data_set:
        out_f = get_sampled_fasta(f)
        fasta_to_token_fasta(ibltseq_exec, toset_exec, "syncmers", f, out_f, args.k, args.z, 1, args.canonical) #sampling rate = 1 because we use minHash to do the heavy lifting
    sampled_query_set = list(map(get_sampled_fasta, query_set))
    sampled_reference_set = list(map(get_sampled_fasta, reference_set))
    mash_sketch(sampled_reference_set, syncmers_mash_ref, args.k, mh_s_param, args.seed, args.canonical)
    mash_dist(syncmers_mash_ref+".msh", sampled_query_set, args.k, mh_s_param, args.seed, tmp_syncmers_mash_file)
    configure4sequences(ibltseq_exec, "syncmers", args.k, args.z, False, False)
    header_names = list()
    for srate in args.rates:
        assert srate >= 0
        header_names.append("syncj r={}".format(srate))
        for f in data_set:
            ibf_sketch_fastx(ibltseq_exec, toset_exec, "syncmers", f, get_ibfon("syncmers", srate, sync_ibf_dir, f), args.k, args.z, ibf_n_param, 3, 0, args.seed, srate, 0, args.canonical)
    # compute similarity scores between Q and R using [exact Jaccard, minHash, syncmers + IBF]
    with open(tmp_other_file, "w") as ohandle:
        ohandle.write("reference,query,size,s,n,exj,{}\n".format(','.join(header_names)))
        for query in query_set:
            for reference in reference_set:
                _, _, _, _, refkmc = get_kmc_paths(args.k, reference, kmc_dataset_dir) #retrieve reference file name
                _, _, _, _, qrykmc = get_kmc_paths(args.k, query, kmc_dataset_dir) #retrieve query file name
                exj, _, _ = kmc_jaccard(cws_exec, qrykmc, refkmc) # ignore weighted jaccard and true symmetric difference size
                sync_ibfjs = list()
                for srate in args.rates:
                    sync_ibfj = compute_pairwise_jaccard_from_ibfs(ibltseq_exec, get_ibfon("syncmers", srate, sync_ibf_dir, query), get_ibfon("syncmers", srate, sync_ibf_dir, reference))
                    sync_ibfjs.append(sync_ibfj)
                ohandle.write("{},{},{},{},{},{},{}\n".format(reference, query, args.ssize, mh_s_param, ibf_n_param, exj, ','.join(map(str, sync_ibfjs)))) #removed smpl_ibfj from csv file
    mash_df = pd.read_csv(tmp_mash_file, sep=',')
    syncmers_mash_df = pd.read_csv(tmp_syncmers_mash_file, sep=',')
    # grouped_df = pd.read_csv(tmp_grouped_syncmers_file, sep=',')
    other_df = pd.read_csv(tmp_other_file, sep=',')
    syncmers_mash_df.rename(columns = {"mhj":"syncmhj"}, inplace = True)
    syncmers_mash_df["reference"] = syncmers_mash_df["reference"].apply(lambda x: x.replace(x, os.path.join(args.ifolder, os.path.basename(x))))
    syncmers_mash_df["query"] = syncmers_mash_df["query"].apply(lambda x: x.replace(x, os.path.join(args.ifolder, os.path.basename(x))))
    joined = pd.merge(mash_df, syncmers_mash_df, on=["reference", "query"], how="inner")
    # joined = pd.merge(joined, grouped_df, on=["reference", "query"], how="inner")
    joined = pd.merge(joined, other_df, on=["reference", "query"], how="inner")# [Post-processing] join results
    joined = joined[["reference", "query", "size", "s", "n", "exj", "mhj", "syncmhj"] + header_names]#, "extendedj", "grouped size"]]
    joined.to_csv(args.o, sep=',', header=True, index=False)

def experiment4_main(args):
    '''Store multiple minHash sketches inside IBFs in order to quickly retrieve sketches associated to interesting sequences for set reconciliation'''
    aldiff_path = os.path.dirname(args.__exepath)
    ibltseq_exec = os.path.join(aldiff_path, main_exec)
    configure4hashes(ibltseq_exec, args.hashwidth, args.s)
    if os.path.isdir(args.alice):
        args.alice = get_fastx_in_folder(args.alice)
    if os.path.isdir(args.bob):
        args.bob = get_fastx_in_folder(args.bob)
    ibf_size = ibf_build_minhash_collection(ibltseq_exec, args.hashwidth, args.s, args.n, args.alice, args.bob)
    print("ibf sketch size is ", ibf_size, " Bytes")

def experiment5_main(args):
    ''''''
    assert args.l >= 0
    assert 0 < args.k
    assert 0 <= args.z < args.k
    assert args.trials > 0
    random.seed(args.seed)
    totj, totn, totm, tots = 0, 0, 0, 0
    sq_errn, sq_errm, sq_errs = 0, 0, 0
    for _ in range(args.trials):
        hash_seed = random.randint(0, 2**64-1)
        j, n, m, s = experiment5(args.__exepath, hash_seed, args.k, args.z)#, seq1, seq2)
        totj += j
        totn += n
        totm += m
        tots += s
        errn = n-j
        errm = m-j
        errs = s-j
        sq_errn += errn**2
        sq_errm += errm**2
        sq_errs += errs**2
    if not args.o: 
        oh = sys.stdout
        header=True
    else: 
        header=isempty(args.o)
        oh = open(args.o, "a") 
    if header: oh.write("L,k,z,mutp,indelf,extp,trials,exact jaccard,sampled jaccard,minimizers jaccard,syncmers jaccard,MSE sampled,MSE minimizers,MSE syncmers,var sampled,var minimizers,var syncmers\n")
    varn = (sq_errn - errn**2 / args.trials) / args.trials
    varm = (sq_errm - errm**2 / args.trials) / args.trials
    vars = (sq_errs - errs**2 / args.trials) / args.trials
    totj = totj/args.trials
    totn = totn/args.trials
    totm = totm/args.trials
    tots = tots/args.trials
    sq_errn = sq_errn/args.trials
    sq_errm = sq_errm/args.trials
    sq_errs = sq_errs/args.trials
    oh.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(args.l, args.k, args.z, args.mutationp, args.indelf, args.extensionp, args.trials, totj, totn, totm, tots, sq_errn, sq_errm, sq_errs,varn,varm,vars)) #CSV header
    oh.close()

def experiment5alt_main(args):
    '''Same thing as experiment5_main except for the output csv'''
    assert args.l >= 0
    assert 0 < args.k
    assert 0 <= args.z < args.k
    assert args.trials > 0
    random.seed(args.seed)
    if not args.o: 
        oh = sys.stdout
        header=True
    else: 
        header=isempty(args.o)
        oh = open(args.o, "a") 
    if header: oh.write("L,k,z,mutp,indelf,extp,trials,exact jaccard,sampled jaccard,minimizers jaccard,syncmers jaccard\n")
    for _ in range(args.trials):
        hash_seed = random.randint(0, 2**64-1)
        truej, sampledj, mmj, syncj = experiment5(args.__exepath, hash_seed, args.k, args.z)#, seq1, seq2)
        oh.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(args.l, args.k, args.z, args.mutationp, args.indelf, args.extensionp, args.trials, truej, sampledj, mmj, syncj)) #CSV header
    oh.close()

def get_kmers(s: str, k: int):
    if len(s) < k: yield s
    for i in range(len(s)-k+1):
        yield s[i:(i+k)]

def get_canonical_kmers(s: str, k: int):
    if len(s) < k: yield s
    for i in range(len(s)-k+1):
        yield kmer.canonical(s[i:(i+k)])

def parse_aggregated_syncmers(multiline_output: str, k: int, canonical: bool) -> set:
    parsed = set()
    for agg in multiline_output:
        if canonical: kmers = get_canonical_kmers(agg, k)
        else: kmers = get_kmers(agg, k)
        parsed |= set(kmers)
    return parsed

def parse_aggregated_symmetric_difference(multiline_output: str, k: int, canonical: bool) -> tuple[set]:
    iset = set()
    jset = set()
    for line in multiline_output:
        source,agg = line.split(',')
        if canonical: kmers = get_canonical_kmers(agg, k)
        else: kmers = get_kmers(agg, k)
        if source == 'i': iset |= set(kmers)
        elif source == 'j': jset |= set(kmers)
        else: raise ValueError("unrecognized input")
    return iset, jset

def compute_exact_symmetric_difference_from_extended_ibfs(executable: str, sk1: str, sk2: str, k: int, canonical: bool, tmpdiff: str) -> int:
    diff_command = [executable, "diff", "-i", sk1, "-j", sk2, "-o", tmpdiff]
    out = subprocess.run(diff_command, stdout=subprocess.DEVNULL, stderr=sys.stderr)
    if out.returncode != 0:
        sys.stderr.write("Error while computing sketch difference between:\n{}\n{}\n".format(sk1, sk2))
        return None
    list_command = [executable, "list", "-i", tmpdiff]
    out = subprocess.run(list_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if out.returncode != 0:
        # sys.stderr.write("Error while peeling the difference sketch of\n{}\n{}\n".format(sk1, sk2))
        return None
    si, sj = parse_aggregated_symmetric_difference(out.stdout.decode("utf-8").splitlines(), k, canonical)
    # count_command = [executable, "count", "-i"]
    # out = subprocess.run(count_command + [sk1], stdout=subprocess.PIPE, stderr=sys.stderr, check=True)
    # L0i = int(out.stdout.decode("utf-8"))
    # out = subprocess.run(count_command + [sk2], stdout=subprocess.PIPE, stderr=sys.stderr, check=True)
    # L0j = int(out.stdout.decode("utf-8"))
    # if sk1 != sk2:
    #     with open("set1.txt", "w") as dbh:
    #         dbh.write('\n'.join(si))
    #     with open("set2.txt", "w") as dbh:
    #         dbh.write('\n'.join(sj))
    return len(si.symmetric_difference(sj))
    # kmbl = number_of_kmers_in_extended_syncmers = (2*k - 2*z + 1) 
    # true_L0i = L0i * kmbl
    # true_L0j = L0j * kmbl
    # number_of_duplicated_kmer, den = kmer.jaccard(si, sj)
    # symm_diff_len = den - number_of_duplicated_kmer
    # union_size = true_L0i + true_L0j - 
    # return num/den, num, den

def experiment6_main(args):
    aldiff_path = os.path.dirname(args.__exepath)
    ibltseq_exec = os.path.join(aldiff_path, main_exec)
    # toset_exec = os.path.join(aldiff_path, toset_script)
    # cws_exec = os.path.join(aldiff_path, kmc_cws_exec)
    random.seed(args.seed)
    seq1 = ''.join(random.choice('ACTG') for _ in range(args.l))
    seq2, _, _, _, _ = mutate.mutate(seq1, args.mutationp, 0, 0)
    grouped_syncmers_command = [ibltseq_exec, "syncmers", "-k", str(args.k), "-m", str(args.z), "-g", str(args.k - args.z)]
    group_out = subprocess.run(grouped_syncmers_command, text=True, input=">0\n{}".format(seq1), stdout=subprocess.PIPE, stderr=sys.stderr, check=True)
    iset = parse_aggregated_syncmers(group_out.stdout.splitlines(), args.k, False)
    if args.d:
        with open(os.path.join(args.d, "grouped_isyncmers.txt"), "w") as oh:
            oh.write("{}".format(group_out.stdout))
    group_out = subprocess.run(grouped_syncmers_command, text=True, input=">0\n{}".format(seq2), stdout=subprocess.PIPE, stderr=sys.stderr, check=True)
    jset = parse_aggregated_syncmers(group_out.stdout.splitlines(), args.k, False)
    if args.d:
        with open(os.path.join(args.d, "grouped_jsyncmers.txt"), "w") as oh:
            oh.write("{}".format(group_out.stdout))
    extended_symm_diff = set(kmer.diff(iset, jset, True))
    if args.d:
        with open(os.path.join(args.d, "seq1.txt"), "w") as oh:
            oh.write("{}\n".format(seq1))
        with open(os.path.join(args.d, "seq2.txt"), "w") as oh:
            oh.write("{}\n".format(seq2))
        with open(os.path.join(args.d, "grouped_iset.txt"), "w") as oh:
            for sync in sorted(list(iset)):
                oh.write("{}\n".format(sync))
        with open(os.path.join(args.d, "grouped_jset.txt"), "w") as oh:
            for sync in sorted(list(jset)):
                oh.write("{}\n".format(sync))
        with open(os.path.join(args.d, "ext_diff.txt"), "w") as oh:
            for sync in sorted(list(extended_symm_diff)):
                oh.write("{}\n".format(sync))
    iset = kmer.set(seq1, args.k, False, None)
    jset = kmer.set(seq2, args.k, False, None)
    exact_symm_diff = set(kmer.diff(iset, jset, True))
    contained = exact_symm_diff.issubset(extended_symm_diff)
    if args.d:
        with open(os.path.join(args.d, "iset.txt"), "w") as oh:
            for sync in sorted(list(iset)):
                oh.write("{}\n".format(sync))
        with open(os.path.join(args.d, "jset.txt"), "w") as oh:
            for sync in sorted(list(jset)):
                oh.write("{}\n".format(sync))
        with open(os.path.join(args.d, "exact_diff.txt"), "w") as oh:
            for sync in sorted(list(exact_symm_diff)):
                oh.write("{}\n".format(sync))
    if not contained: 
        print("Not ok")
        print(extended_symm_diff - exact_symm_diff)
        print(exact_symm_diff - extended_symm_diff)
    else: print("OK")

def experiment7_main(args):
    '''I know, code duplication from experiment 3 is bad but I got no time to refactor.
    '''
    assert 0 < args.k <= 32
    assert 0 <= args.z <= args.k
    assert (args.qnumber and args.qnumber >= 0) or (args.qfraction and 0 <= args.qfraction <= 1)
    assert os.path.isdir(args.ifolder)
    assert args.epsilon >= 0
    aldiff_path = os.path.dirname(args.__exepath)
    ibltseq_exec = os.path.join(aldiff_path, main_exec)
    toset_exec = os.path.join(aldiff_path, toset_script)
    cws_exec = os.path.join(aldiff_path, kmc_cws_exec)

    random.seed(args.seed)
    kmc_dir = os.path.join(args.wfolder, "kmc")
    _, dataset_basename = get_path_leaf(args.ifolder)
    kmc_dataset_dir = os.path.join(kmc_dir, dataset_basename) # used later to retrieve kmc databases
    kmc_tmp_dir = os.path.join(args.wfolder, "kmc_tmp")
    sync_ibf_dir = os.path.join(args.wfolder, "syncmer_ibf")
    os.makedirs(kmc_dir, exist_ok=True)
    os.makedirs(kmc_tmp_dir, exist_ok=True)
    os.makedirs(sync_ibf_dir, exist_ok=True)

    payload_size = int((args.k * 2 + 8 - 1) / 8)
    counter_size = 8
    bucket_size = payload_size + counter_size
    repetitions = 3
    peelability_oversize = 1.222 + args.epsilon
    def get_sketch_size(actual_n: int):
        return int(round(bucket_size * peelability_oversize * actual_n + bucket_size * repetitions))
    
    data_set = get_fastx_in_folder(args.ifolder)
    if args.qnumber and not args.qfraction: query_size = args.qnumber
    elif not args.qnumber and args.qfraction: query_size = max(1, int(len(data_set) * args.qfraction))
    else: raise RuntimeError("Something went wrong with the mutually-exclusive options 'qnumber' and 'qfraction'")
    query_set = sorted(random.sample(data_set, query_size))
    if (query_size == len(data_set)): reference_set = data_set
    else: reference_set = sorted(list(set(data_set) - set(query_set)))

    _ = kmc_count_folder(args.ifolder, kmc_dir, [args.k], kmc_tmp_dir, 8, False, args.canonical)# Count files using KMC
    #find maximum n and build kmc dataframe
    def insert_kmc_row(df, reference: str, query: str, exact_jaccard: float):
        df.loc[len(df)] = [reference, query, exact_jaccard]
    mn = 0
    kmc_df = pd.DataFrame(columns=["reference", "query", "exj"])
    for query in query_set:
        for reference in reference_set:
            _, _, _, _, refkmc = get_kmc_paths(args.k, reference, kmc_dataset_dir) #retrieve reference file name
            _, _, _, _, qrykmc = get_kmc_paths(args.k, query, kmc_dataset_dir) #retrieve query file name
            exj, _, true_n = kmc_jaccard(cws_exec, qrykmc, refkmc) # ignore weighted jaccard
            insert_kmc_row(kmc_df, reference, query, exj)
            true_n = int(true_n)
            mn = max(mn, true_n)
    configure4sequences(ibltseq_exec, "syncmers", args.k, args.z, False, False)
    def get_sampled_n(maximum_n: int, k: int, z: int, sampling_rate: int):
        return int(maximum_n / ((k-z+1)/2 * sampling_rate))
    for srate in args.rates:
        assert srate >= 0
        for f in data_set:
            sketch_name = get_ibfon("syncmers", srate, sync_ibf_dir, f)
            # if not os.path.isfile(sketch_name):
            sampled_n = get_sampled_n(mn, args.k, args.z, srate)
            if sampled_n == 0: raise RuntimeError("IBF sketches of size 0 given the current parameters")
            ibf_sketch_fastx(ibltseq_exec, toset_exec, "syncmers", f, sketch_name, args.k, args.z, sampled_n, 3, args.epsilon, args.seed, srate, 0, args.canonical)

    # compute exact Jaccard, and estimations using sampled syncmers for different n parameters
    def insert_sampled_row(df, rate: int, reference: str, query: str, size: int, syncmer_jaccard: float):
        df.loc[len(df)] = [rate, reference, query, size, syncmer_jaccard]
    sync_df = pd.DataFrame(columns=["sampling rate","reference","query","size","syncj"])
    for srate in args.rates:
        sampled_n = get_sampled_n(mn, args.k, args.z, srate)
        sketch_size = get_sketch_size(sampled_n)
        for query in query_set:
            ibf_query_name = get_ibfon("syncmers", srate, sync_ibf_dir, query)
            for reference in reference_set:
                ibf_ref_name = get_ibfon("syncmers", srate, sync_ibf_dir, reference)
                syncj = compute_pairwise_jaccard_from_ibfs(ibltseq_exec, ibf_query_name, ibf_ref_name)
                insert_sampled_row(sync_df, srate, reference, query, sketch_size, syncj)
                #ohandle.write("{},{},{},{},{},{},{}\n".format(srate, reference, query, sketch_size, syncj)) #removed smpl_ibfj from csv file
    joined = pd.merge(kmc_df, sync_df, on=["reference", "query"], how="inner")# [Post-processing] join results
    if os.path.isfile(args.o) and os.stat(args.o).st_size != 0:
        source_df = pd.read_csv(args.o, sep=',')
        joined = source_df.append(joined)
    joined.to_csv(args.o, sep=',', header=True, index=False)

def fast2set(fastx_file: str, k: int, canonical: bool) -> set:
    table = set()
    with open(fastx_file, "r") as fi:
        for _, seq, _ in fastx.read(fi):
            kmer.set(seq, k, canonical, table)
    return table

def ibf_sketch_fastx_alt(executable: str, toset: str, algorithm: str, input_file: str, output_file: str, k: int, z: int, n: int, r: int, epsilon: float, seed: int, extra_len: int, canonical: bool):
    """Alternative version for Ex 8"""
    fragmentation_command = [executable]
    if algorithm == "syncmers": fragmentation_command.extend(["syncmers", "-m", str(z), "-g", str(extra_len)])
    elif algorithm == "minimizers": fragmentation_command.extend(["minimizers", "-m", str(z)])
    elif algorithm == "kmers": fragmentation_command.extend(["kmers"])
    else: sys.stderr.write("This should never happen\n") #This should never happen because configure_length raise error if algorithm unrecognized
    fragmentation_command.extend(["-i", input_file, "-k", str(k)])
    toset_command = ["python3", toset]
    build_command = [executable, "build", "-o", output_file, "-n", str(n), "-r", str(r), "-e", str(epsilon), "-s", str(seed)]
    fragmentation_out = subprocess.Popen(fragmentation_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if canonical: toset_command.append("-c")
    toset_out = subprocess.Popen(toset_command, stdin=fragmentation_out.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    build_out = subprocess.Popen(build_command, stdin=toset_out.stdout, stdout=sys.stdout, stderr=sys.stderr)
    build_out.communicate()
    toset_out.communicate()#just to be sure
    if build_out.returncode != 0:
        sys.stderr.write("building from fastx file exited with error\n")

def experiment8_main(args):
    assert 0 < args.k <= 32
    assert 0 <= args.z <= args.k
    # assert (args.qnumber and args.qnumber >= 0) or (args.qfraction and 0 <= args.qfraction <= 1)
    assert args.ssize >= 0
    assert os.path.isdir(args.ifolder)
    aldiff_path = os.path.dirname(args.__exepath)
    ibltseq_exec = os.path.join(aldiff_path, main_exec)
    toset_exec = os.path.join(aldiff_path, toset_script)
    cws_exec = os.path.join(aldiff_path, kmc_cws_exec)

    random.seed(args.seed)
    tmp_diff_ibf = os.path.join(args.wfolder, get_random_name("ibf.bin"))
    extended_ibf_dir = os.path.join(args.wfolder, "extended_syncmers_ibf")
    kmc_dir = os.path.join(args.wfolder, "kmc")
    _, dataset_basename = get_path_leaf(args.ifolder)
    kmc_dataset_dir = os.path.join(kmc_dir, dataset_basename) # used later to retrieve kmc databases
    kmc_tmp_dir = os.path.join(args.wfolder, "kmc_tmp")
    extended_ibf_dir = os.path.join(args.wfolder, "extended_syncmers_ibf")
    os.makedirs(kmc_dir, exist_ok=True)
    os.makedirs(kmc_tmp_dir, exist_ok=True)
    os.makedirs(extended_ibf_dir, exist_ok=True)

    klen = args.k + 1 * (args.k - args.z)
    payload_size = int((klen * 2 + 8 - 1) / 8)
    counter_size = 8
    length_cell_size = 0 #1 Byte -> up to lenngth 256. Since we limit k to be <=32 everything is ok
    bucket_size = payload_size + counter_size + length_cell_size
    repetitions = 3
    epsilon = 0
    peelability_oversize = 1.222 + epsilon
    ibf_n_param = math.ceil((args.ssize - bucket_size * repetitions) / (bucket_size * peelability_oversize)) # maximum number of differences trackable by the given sketch size
    if ibf_n_param == 0: raise RuntimeError("IBF sketches of size 0 given the current parameters")
    data_set = get_fastx_in_folder(args.ifolder)
    # Count files using KMC
    _ = kmc_count_folder(args.ifolder, kmc_dir, [args.k], kmc_tmp_dir, 8, False, args.canonical)
    reference_set = query_set = data_set
    configure4sequences(ibltseq_exec, "syncmers", klen, args.z,False, False)
    for f in data_set:
        ibf_sketch_fastx_alt(ibltseq_exec, toset_exec, "syncmers", f, get_ibfon("syncmers", 1, extended_ibf_dir, f), args.k, args.z, ibf_n_param, 3, 0, args.seed, args.k - args.z, args.xanonical)
    def insert_extended_row(df, reference: str, query: str, size: int, n: int, true_symmd, est_symmd):
        df.loc[len(df)] = [reference, query, size, n, true_symmd, est_symmd]
    sync_df = pd.DataFrame(columns=["reference","query","size","n","true symmetric size","estimated symmetric size"])
    for query in query_set:
        for reference in reference_set:
            _, _, _, _, refkmc = get_kmc_paths(args.k, reference, kmc_dataset_dir) #retrieve reference file name
            _, _, _, _, qrykmc = get_kmc_paths(args.k, query, kmc_dataset_dir) #retrieve query file name
            exj, _, true_difference_size = kmc_jaccard(cws_exec, qrykmc, refkmc) # ignore weighted jaccard
            estimated_difference_size = compute_exact_symmetric_difference_from_extended_ibfs(ibltseq_exec, get_ibfon("syncmers", 1, extended_ibf_dir, query), get_ibfon("syncmers", 1, extended_ibf_dir, reference), args.k, args.canonical, tmp_diff_ibf)
            insert_extended_row(sync_df, reference, query, args.ssize, ibf_n_param, int(true_difference_size), estimated_difference_size)
    sync_df.to_csv(args.o, sep=',', header=True, index=False)

def main(args):
    if args.command == "sketch": return sketch_main(args)
    elif args.command == "multi": return multi_main(args)
    elif args.command == "jaccard": return jaccard_main(args)
    elif args.command == "rndset": return rndset_main(args)
    elif args.command == "ex1": return experiment1_main(args)
    elif args.command == "ex2": return experiment2_main(args)
    elif args.command == "ex3": return experiment3_main(args)
    elif args.command == "ex4": return experiment4_main(args)
    elif args.command == "ex5": return experiment5_main(args)
    elif args.command == "ex5alt": return experiment5alt_main(args)
    elif args.command == "ex6": return experiment6_main(args)
    elif args.command == "ex7": return experiment7_main(args)
    elif args.command == "ex8": return experiment8_main(args)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_sketch = subparsers.add_parser("sketch", help="Construct IBF from input set of sequences")
    parser_sketch.add_argument("-i", help="input file", type=str, required=True)
    parser_sketch.add_argument("-o", help="output IBF", type=str, required=True)
    parser_sketch.add_argument("-l", help="length of strings in the input set", type=int, required=True)
    parser_sketch.add_argument("-n", help="expected number of differences between two IBFs", type=int, required=True)
    parser_sketch.add_argument("-r", "--repetitions", help="number of hash functions", type=int, default=3)
    parser_sketch.add_argument("-e", "--epsilon", help="over-dimensioning factor of the IBF", type=float, default=0)
    parser_sketch.add_argument("-s", "--seed", help="random seed", type=int, default=42)

    parser_multi = subparsers.add_parser("multi", help="Construct one IBF for each record in a multi-fasta file")
    parser_multi.add_argument("-i", help="input file", type=str, required=True)
    parser_multi.add_argument("-o", help="output folder", type=str, required=True)
    parser_multi.add_argument("-a", "--algorithm", help="algorithm to compute the set (minimizers, syncmers) [syncmers]", type=str, default="syncmers", choices=["minimizers", "syncmers"])
    parser_multi.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_multi.add_argument("-m", help="minimizer length:\n\tif 'minimizers' method is selected m is used to group k-mers together\n\tif 'syncmers' is selected, m is used to define syncmers", type=int, required=True)
    parser_multi.add_argument("-n", help="expected number of differences between two IBFs", type=int, required=True)
    parser_multi.add_argument("-f", help="fragment each sequence at minimizer positions (option -a becomes irrelevant)", action="store_true")
    parser_multi.add_argument("-r", "--repetitions", help="number of hash functions", type=int, default=3)
    parser_multi.add_argument("-e", "--epsilon", help="over-dimensioning factor of the IBF", type=float, default=0)
    parser_multi.add_argument("-s", "--seed", help="random seed", type=int, default=42)

    parser_jaccard = subparsers.add_parser("jaccard", help="Compute pairwise jaccard and containments between multiple sketches")
    parser_jaccard.add_argument("-i", help="list of sketches or folder containing them", type=str, nargs='+', required=True)
    parser_jaccard.add_argument("-o", help="output csv file to append to", type=str, required=True)

    parser_rndset = subparsers.add_parser("rndset", help="Generate random dataset. Each sequence is the mutated version of a random founder")
    parser_rndset.add_argument("-o", help="output prefix (a path in the form of <output folder>/<name prefix>)", type=str, required=True)
    parser_rndset.add_argument("-n", help="number of genomes", type=int, required=True)
    parser_rndset.add_argument("-l", help="length of simulated genomes", type=int, required=True)
    parser_rndset.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    parser_rndset.add_argument("-d", "--indelf", help="fraction of mutations that are indels [0]", type=float, default=0)
    parser_rndset.add_argument("-x", "--extensionp", help="probability of extending an indel [0]", type=float, default=0)
    parser_rndset.add_argument("-s", "--seed", help="main random seed used for sequence generation [42]", type=int, default=42)

    parser_ex1 = subparsers.add_parser("ex1", help="Experiment 1: compare IBFs built using minimizers and syncmers on random uniform sequences")
    parser_ex1.add_argument("-o", help="output csv file to append to [stdout]", type=str)
    parser_ex1.add_argument("-n", help="number of genomes", type=int, required=True)
    parser_ex1.add_argument("-l", help="length of simulated genomes", type=int, required=True)
    parser_ex1.add_argument("-k", help="k-mer (syncmer) length", type=int, required=True)
    parser_ex1.add_argument("-z", help="minimizer length used for grouping or to define syncmers", type=int, required=True)
    parser_ex1.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    parser_ex1.add_argument("-d", "--indelf", help="fraction of mutations that are indels [0]", type=float, default=0)
    parser_ex1.add_argument("-x", "--extensionp", help="probability of extending an indel [0]", type=float, default=0)
    parser_ex1.add_argument("-t", "--mtrials", help="maximum n for each IBF", type=int, default=math.inf)
    parser_ex1.add_argument("-s", "--seed", help="main random seed used for sequence generation [42]", type=int, default=42)

    parser_ex2 = subparsers.add_parser("ex2", help="Experiment 2: measure IBF performances")
    parser_ex2.add_argument("-o", help="output csv file to append to [stdout]", type=str)
    parser_ex2.add_argument("-w", "--wfolder", help="folder used to store temporary files", type=str, required=True)
    parser_ex2.add_argument("-a", "--algorithm", help="algorithm used to generate input sets. (minimizers, syncmers) [syncmers]", type=str, default="syncmers", choices=["minimizers", "syncmers"])
    parser_ex2.add_argument("-n", help="number of genomes", type=int, required=True)
    parser_ex2.add_argument("-l", help="length of simulated genomes", type=int, required=True)
    parser_ex2.add_argument("-k", help="k-mer (syncmer) length", type=int, required=True)
    parser_ex2.add_argument("-z", help="minimizer length used for grouping or to define syncmers", type=int, required=True)
    parser_ex2.add_argument("-e", help="IBF epsilon parameter [0]", type=float, default=0)
    parser_ex2.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    parser_ex2.add_argument("-d", "--indelf", help="fraction of mutations that are indels [0]", type=float, default=0)
    parser_ex2.add_argument("-x", "--extensionp", help="probability of extending an indel [0]", type=float, default=0)
    parser_ex2.add_argument("-s", "--seed", help="main random seed used for sequence generation [42]", type=int, default=42)

    parser_ex3 = subparsers.add_parser("ex3", help="Experiment 3: compare minHash (MASH) and syncmers + IBF")
    parser_ex3.add_argument("-i", "--ifolder", help="folder containing fasta files to be sketched", type=str, required=True)
    parser_ex3.add_argument("-o", help="output csv file", type=str, required=True)
    parser_ex3.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_ex3.add_argument("-z", help="syncmer selection parameter for the syncmers+IBF case", type=int, required=True)
    parser_ex3.add_argument("-m", "--ssize", help="sketch size in Bytes for both minHashes and IBFs", type=int, required=True)
    query_option_parser = parser_ex3.add_mutually_exclusive_group(required = True)
    query_option_parser.add_argument("-q", "--qnumber", help="number of randomly selected sequences to be used as queries (disables -f option). If 0 use all sequences.", type=int)
    query_option_parser.add_argument("-f", "--qfraction", help="fraction of input sequences to be considered as queries (disables -q option). If 0 use all sequences.", type=float)
    parser_ex3.add_argument("-r", "--rates", help="syncmers's sampling rates [1]", type=int, nargs='+', default=[1])
    parser_ex3.add_argument("-w", "--wfolder", help="folder used to create temporary files", type=str, required=True)
    parser_ex3.add_argument("-c", "--canonical", help="activate canonical k-mers", action="store_true")
    parser_ex3.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)

    parser_ex4 = subparsers.add_parser("ex4", help="Experiment 4: combine minHash sketches with IBF to produce small representations of sequence databases to be used for remote database reconciliation")
    parser_ex4.add_argument("-A", "--Alice", help="first collection of sequences", type=str, required=True)
    parser_ex4.add_argument("-B", "--Bob", help="second database of sequences", type=str, required=True)
    parser_ex4.add_argument("-b", "--bitwidth", help="hash width [64]", type=int, default=64)
    parser_ex4.add_argument("-s", help="number of hashes in each minHash sketch", type=int, required=True)
    parser_ex4.add_argument("-n", help="maximum number of dissimilar sketches to track", type=int, required=True)
    parser_ex4.add_argument("-S", "--seed", help="random seed [42]", type=int, default=42)

    parser_ex5 = subparsers.add_parser("ex5", help="Experiment 5: test if syncmers are truly un-biased")
    parser_ex5.add_argument("-o", help="output csv", type=str, required=False)
    parser_ex5.add_argument("-l", help="length of the simulated genomes", type=int, required=True)
    parser_ex5.add_argument("-k", help="k-mer (syncmer) length", type=int, required=True)
    parser_ex5.add_argument("-z", help="minimizer length used for grouping or to define syncmers", type=int,  required=True)
    parser_ex5.add_argument("-t", "--trials", help="number of trials (different random sequences)", type=int, required=True)
    parser_ex5.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    parser_ex5.add_argument("-d", "--indelf", help="fraction of mutations that are indels [0]", type=float, default=0)
    parser_ex5.add_argument("-x", "--extensionp", help="probability of extending an indel [0]", type=float, default=0)
    parser_ex5.add_argument("-s", "--seed", help="main random seed used for sequence generation [42]", type=int, default=42)

    parser_ex5alt = subparsers.add_parser("ex5alt", help="Experiment 5: test if syncmers are truly un-biased, alternative version outputting all measures")
    parser_ex5alt.add_argument("-o", help="output csv", type=str, required=False)
    parser_ex5alt.add_argument("-l", help="length of the simulated genomes", type=int, required=True)
    parser_ex5alt.add_argument("-k", help="k-mer (syncmer) length", type=int, required=True)
    parser_ex5alt.add_argument("-z", help="minimizer length used for grouping or to define syncmers", type=int,  required=True)
    parser_ex5alt.add_argument("-t", "--trials", help="number of trials (different random sequences)", type=int, required=True)
    parser_ex5alt.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    parser_ex5alt.add_argument("-d", "--indelf", help="fraction of mutations that are indels [0]", type=float, default=0)
    parser_ex5alt.add_argument("-x", "--extensionp", help="probability of extending an indel [0]", type=float, default=0)
    parser_ex5alt.add_argument("-s", "--seed", help="main random seed used for sequence generation [42]", type=int, default=42)
    
    parser_ex6 = subparsers.add_parser("ex6", help="Extended syncmers on simulated sequences: check if true supersets and compute Jaccard")
    parser_ex6.add_argument("-d", help="output prefix for debugging", type=str, required=False)
    parser_ex6.add_argument("-l", help="length of simulated genomes", type=int, required=True)
    parser_ex6.add_argument("-k", help="k-mer (syncmer) length", type=int, required=True)
    parser_ex6.add_argument("-z", help="syncmer selection parameter for syncmers+IBF case", type=int,  required=True)
    parser_ex6.add_argument("-m", "--mutationp", help="probability of mutation", type=float, required=True)
    # parser_ex6.add_argument("-g", help="aggregate g neighboring k-mers to each syncmer", type=int, required=True)
    parser_ex6.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)

    parser_ex7 = subparsers.add_parser("ex7", help="Compare sampled syncmers + IBLTs space")
    parser_ex7.add_argument("-i", "--ifolder", help="folder containing fasta files to be sketched", type=str, required=True)
    parser_ex7.add_argument("-o", help="output csv file to append to", type=str, required=True)
    parser_ex7.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_ex7.add_argument("-z", help="syncmer selection parameter for the syncmers+IBF case", type=int, required=True)
    parser_ex7.add_argument("-e", "--epsilon", help="IBF epsilon parameter [0]", type=float, default=0)
    query_option_parser = parser_ex7.add_mutually_exclusive_group(required = True)
    query_option_parser.add_argument("-q", "--qnumber", help="number of randomly selected sequences to be used as queries (disables -f option). If 0 use all sequences.", type=int)
    query_option_parser.add_argument("-f", "--qfraction", help="fraction of input sequences to be considered as queries (disables -q option). If 0 use all sequences.", type=float)
    parser_ex7.add_argument("-r", "--rates", help="sampling rates of the combination of syncmers and sampling [1]", type=int, nargs='+', default=[1])
    parser_ex7.add_argument("-w", "--wfolder", help="folder used to create temporary files", type=str, required=True)
    parser_ex7.add_argument("-c", "--canonical", help="activate canonical k-mers", action="store_true")
    parser_ex7.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)

    parser_ex8 = subparsers.add_parser("ex8", help="Extended syncmers on real datasets: try to estimate jaccard")
    parser_ex8.add_argument("-i", "--ifolder", help="folder containing fasta files to be sketched", type=str, required=True)
    parser_ex8.add_argument("-o", help="output csv file", type=str, required=True)
    parser_ex8.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_ex8.add_argument("-z", help="syncmer selection parameter for the syncmers+IBF case", type=int, required=True)
    parser_ex8.add_argument("-m", "--ssize", help="sketch size in Bytes for both minHashes and IBFs", type=int, required=True)
    parser_ex8.add_argument("-c", "--canonical", help="activate canonical disaggregated k-mers", action="store_true")
    parser_ex8.add_argument("-x", "--xanonical", help="store canonical extended segments inside IBFs", action="store_true")
    parser_ex8.add_argument("-w", "--wfolder", help="folder used to create temporary files", type=str, required=True)
    parser_ex8.add_argument("-s", "--seed", help="random seed [42]", type=int, default=42)

    return parser

if __name__ == "__main__":
    mypath = os.path.dirname(sys.argv[0])
    abspath = os.path.abspath(mypath)
    parser = parser_init()
    args = parser.parse_args(sys.argv)
    setattr(args, "__exepath", abspath)
    main(args)

def get_data_set(folder: str):
    """use get_fastx_in_folder() instead
    """
    data_set = list()
    for f in os.listdir(folder):
        full_f_path = os.path.join(folder,f)
        if os.path.isfile(full_f_path):
            try:
                get_fastx_format(full_f_path)
                data_set.append(full_f_path)
            except ValueError:
                None
    data_set.sort()

