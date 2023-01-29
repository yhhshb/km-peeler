import math
import pathlib
import subprocess
import multiprocessing

def get_size_from_byte_size(byte_size: int, k: int) -> int:
    assert byte_size >= 0
    assert k >= 0
    return math.ceil(byte_size / (4 if 4**k < 2**32 else 8)) # number of hashes in a minHash sketch (automatic selection between 4 or 8 bytes)

def mash_sketch(reference_set: list[str], output_path: str, k: int, s: int, seed: int, canonical: bool, force: bool = True):
    # use mash sketch of multiple sequences + mash dist of multiple files
    if not force and pathlib.Path(output_path).with_suffix(".msh").exists(): return
    sketch_command = ["mash", "sketch", "-o", output_path, "-k", str(k), "-s", str(s), "-S", str(seed)] #ignore strand
    if not canonical: sketch_command.extend(["-n"])
    sketch_command.extend(reference_set)
    out = subprocess.run(sketch_command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if out.returncode != 0:
        raise RuntimeError("Error, mash sketch command returned with exit code {}\n".format(out.returncode))
        # sys.exit(os.EX_CANTCREAT)

def mash_dist(reference_database, queries: list[str], k: int, s: int, seed: int, tmp_file, force: bool = True):
    if not force and pathlib.Path(tmp_file).exists(): return
    reference_database = pathlib.Path(reference_database)
    dist_command = ["mash", "dist", "-s", str(s), "-S", str(seed), "-p", str(multiprocessing.cpu_count())]
    if reference_database.suffix != ".msh": dist_command.extend(["-k", str(k)])
    dist_command.append(reference_database)
    dist_command.extend(queries)
    out = subprocess.run(dist_command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if out.returncode != 0:
        raise RuntimeError("Error from mash dist command\n")
        # sys.exit(os.EX_CANTCREAT)
    lines = out.stdout.decode("utf-8").split('\n')
    with open(tmp_file, "w") as mhandle:
        mhandle.write("reference,query,mhj\n")
        for line in lines:
            if line:
                rf, qr, _, _, frac = line.split('\t')
                n, d = frac.split('/')
                j = int(n)/int(d)
                mhandle.write("{},{},{}\n".format(rf,qr,j))