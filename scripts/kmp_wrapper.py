import os
import sys
import math
import pathlib
import subprocess
import essentials

ck_table = [0, 0, 0, 1.222, 1.295, 1.425, 1.570, 1.721]

def get_bit_size(n: int, k: int, r: int, epsilon: float):
    bit_counter_size = 2 # 2 bits per counter in the current implementation
    hrc_bit_size = math.ceil((r - 2) * math.log(n, 2) + r)
    bit_payload_size = int((2*k + hrc_bit_size + 7) / 8) * 8
    bucket_size = bit_counter_size + bit_payload_size
    chunk_len = int((ck_table[r] + epsilon) * n / r) + 1
    return r * chunk_len * bucket_size
    
def get_maximum_difference_from_byte_size(byte_size: int, k: int, r: int, epsilon: float):
    assert byte_size >= 0
    assert k >= 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    bit_size = byte_size * 8
    new_guess = int(bit_size / get_bit_size(1, k, r, epsilon))
    guess = [0 for _ in range(3)]
    i = 0
    positive_front = 0
    while(new_guess not in guess and positive_front != new_guess):
        guess[i] = new_guess
        # print(guess)
        size = get_bit_size(guess[i], k, r, epsilon)
        if size > bit_size:
            positive_front = guess[i]
            size_overshoot = (size - bit_size) / bit_size
            # print("-" + str(size_overshoot))
            new_guess = int(guess[i] - guess[i] * size_overshoot)
        else:
            size_undershoot = (bit_size - size) / bit_size
            # print("+" + str(size_undershoot))
            new_guess = int(guess[i] + guess[i] * size_undershoot)
        i = (i + 1) % 3
    return new_guess

def sketch(executable, input_fastx: str, output_sketch: str, n: int, k: int, z: int, x: int, r: int, epsilon: float, seed: int, canonical: bool, tmp_dir, max_ram: int, log_file: str):
    assert executable
    assert input_fastx
    assert output_sketch
    assert 0 < k <= 32
    assert 0 <= z <= k
    assert n > 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    assert seed > 0
    assert tmp_dir
    assert max_ram > 0
    build = [str(executable), "build"]
    input_opt = [str(input_fastx)]
    output_opt = [str(output_sketch)]
    n_opt = ["-n", str(n)]
    k_opt = ["-k", str(k)]
    z_opt = ["-z", str(z)]
    x_opt = ["-x", str(x)]
    r_opt = ["-r", str(r)]
    epsilon_opt = ["-e", str(epsilon)]
    seed_opt = ["-s", str(seed)]
    canon_opt = ["-c"] if canonical else []
    tmp_opt = ["-d", tmp_dir]
    logging_opt = ["--log", log_file] if log_file else []
    command = build + input_opt + output_opt + n_opt + k_opt + z_opt + x_opt + r_opt + epsilon_opt + seed_opt + canon_opt + tmp_opt + logging_opt
    out = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if out.returncode != 0: sys.stderr.write("ibltseq exited with error {}\n".format(out.returncode))

def get_outname(fastx_file) -> str:
    p = pathlib.Path(fastx_file)
    while(p.stem != str(p)): p = pathlib.Path(p.stem)
    return str(p.with_suffix(".iblt.bin"))
    
def sketch_folder(executable, input_folder, output_folder, n: int, k: int, z: int, r: int, epsilon: float, seed: int, canonical: bool, tmp_dir, max_ram:int, force: bool = True):
    assert executable
    assert input_folder
    assert output_folder
    assert 0 < k <= 32
    assert 0 <= z <= k
    assert n > 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    assert seed > 0
    assert tmp_dir
    assert max_ram > 0
    output_folder = pathlib.Path(output_folder)
    fastx_files = essentials.get_fastx_in_folder(input_folder)
    for fname in fastx_files:
        oname = output_folder.joinpath(get_outname(fname))
        if (force or not oname.exists()):
            sketch(executable, fname, oname, n, k, z, 0, r, epsilon, seed, canonical, tmp_dir, max_ram, "")

def pairwise_jaccard(executable, a, b) -> str:
    out = subprocess.run([executable, "diff", str(a), str(b), "-j", "."], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if out.returncode != 0: raise RuntimeError("IBLT jaccard computation failed")
    return out.stdout.decode("utf-8").strip().split(',')[-1]

def check_enhanced_extended_syncmers(executable, a, b, ecc_a, ecc_b, max_ram: int, A, B):
    out = subprocess.run([executable, "correct", a, ecc_a, b, ecc_b, "-m", str(max_ram), "--check", A, B], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    msg = out.stderr.decode("utf-8")
    if out.returncode != 0: 
        # sys.stderr.write(msg)
        raise RuntimeError("Unable to remove false positives")
    return msg.find("OK") == (len(msg) - 3)

def sketch_each_sequence(executable, algo: str, input_fastx, output_folder, k: int, m: str, n: int, r: int, epsilon: float, seed: int):
    #TODO stream temporary fasta to kmp build directly, without writing to disk
    assert input_fastx
    assert output_folder and os.path.isdir(output_folder)
    assert 0 < k <= 32
    assert 0 <= m < k
    assert n > 0
    assert 3 <= r <= 7
    assert 0 <= epsilon <= 1
    executable = pathlib.Path(executable)
    aldiff_path = executable.parent
    build = [executable, "build"]
    k_opt = ["-k", str(k)]
    n_opt = ["-n", str(n)]
    r_opt = ["-r", str(r)]
    epsilon_opt = ["-e", str(epsilon)]
    seed_opt = ["-s", str(seed)]
    partial_build_command = build + k_opt + n_opt + r_opt + epsilon_opt + seed_opt
    input_basename = pathlib.Path(input_fastx).stem
    dummy_fasta_name = essentials.get_random_name("fna")
    dummy_input_name = essentials.get_random_name("txt")
    with open(dummy_fasta_name, "w") as dh:
        with open(input_fastx, "r") as fh:
            for name, seq, _ in fastx.read(fh):
                dh.seek(0)
                dh.write(">{}\n{}".format(name, seq))
                dh.truncate()
                dh.flush()
                fragment_command = "{} {} -i {} -k {} -m {} | python3 {} > {}".format(executable, algo, dummy_fasta_name, k, m, aldiff_path.joinpath(toset_script), dummy_input_name)
                out = subprocess.run(fragment_command, capture_output=False, shell=True)
                if out.returncode != 0: 
                    sys.stderr.write("Unable to compute {} table for record {} (failed with exit code {})\n".format(algo, name, out.returncode))
                else:
                    output_name = "{}_{}.{}".format(input_basename, name, "ibf.bin")
                    output_opt = ["-o", os.path.join(output_folder, output_name)]
                    input_opt = ["-i", dummy_input_name]
                    out = subprocess.run(partial_build_command + input_opt + output_opt, stdout=sys.stdout, stderr=sys.stderr)
                    if out.returncode != 0: sys.stderr.write("ibltseq construction for record {} failed with exit code {}\n", name, out.returncode)