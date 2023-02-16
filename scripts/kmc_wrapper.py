import os
import sys
import pathlib
import subprocess
import essentials

def get_kmc_paths(k: int, input_file, output_path):
    _, filename = essentials.get_path_leaf(input_file)
    filename = filename.split('.')[0]
    kfolder = pathlib.Path(output_path).joinpath("k{}".format(k))
    intermidiate_file = kfolder.joinpath("dummy_{}".format(filename))
    sorted_file = kfolder.joinpath("S_{}".format(filename))
    final_file = kfolder.joinpath(filename)
    return filename, kfolder, intermidiate_file, sorted_file, final_file

def kmc_count(k: int, input_fastx: str, output_path: str, kmc_dummy_folder: str, mmemory: int, unsorted: bool, canonical: bool, force: bool = True) -> str:
    """Call kmc for k-mer counting"""
    _, kfolder, intermidiate_file, sorted_file, final_file = get_kmc_paths(k, input_fastx, output_path)
    final_prefix_file = final_file.with_suffix(".kmc_pre")
    final_suffix_file = final_file.with_suffix(".kmc_suf")
    if (not force) and final_prefix_file.exists() and final_suffix_file.exists(): return final_file # Avoid recomputation
    fmt, _ = essentials.get_fastx_format(input_fastx)
    if fmt == "fasta": fmt = "-fm"
    elif fmt == "fastq": fmt = "-fq"
    else: raise ValueError("Unrecognized file type")

    try: os.makedirs(kfolder)
    except FileExistsError: pass

    kmc_command = ["kmc"]
    if not canonical: kmc_command.extend(["-b"])
    kmc_command.extend(["-k"+str(k), "-ci0", "-cs4294967295", "-cx4294967295", "-m"+str(mmemory), fmt, input_fastx, str(intermidiate_file), str(kmc_dummy_folder)])
    out = subprocess.run(kmc_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) # on macOS terminal write "ulimit -n 2048" before launching
    if out.returncode != 0: raise RuntimeError("KMC failed on file " + input_fastx)
    if (not unsorted):
        out = subprocess.run(["kmc_tools", "transform", intermidiate_file, "sort", sorted_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if out.returncode != 0: raise RuntimeError("Unable to sort KMC output for " + input_fastx)
    
    intermidiate_prefix_file = intermidiate_file.with_suffix(".kmc_pre")
    sorted_prefix_file = sorted_file.with_suffix(".kmc_pre")
    if sorted_prefix_file.exists(): intermidiate_prefix_file.unlink()
    else: intermidiate_prefix_file.rename(sorted_prefix_file)

    intermidiate_suffix_file = intermidiate_file.with_suffix(".kmc_suf")
    sorted_suffix_file = sorted_file.with_suffix(".kmc_suf")
    if sorted_suffix_file.exists(): intermidiate_suffix_file.unlink()
    else: intermidiate_suffix_file.rename(sorted_suffix_file)

    sorted_prefix_file.replace(final_prefix_file)
    sorted_suffix_file.replace(final_suffix_file)
    return final_file

def kmc_count_folder(input_folder, output_folder, ks: list[int], kmc_dummy_folder: str, mmemory: int, unsorted: bool, canonical: bool, force: bool = True) -> list[str]:
    """Build kmc databases for the input folder for all k values.
    
    The output folder will have a folder of the same name as the input with multiple folders named as k<value>.
    This function also creates the temporary folder for kmc hard-disk computations.
    """
    input_dir = pathlib.Path(input_folder)
    output_dir = pathlib.Path(output_folder)
    assert input_dir.is_dir()
    assert type(ks) == type(list())
    assert 0 < mmemory
    fastxs = essentials.get_fastx_in_folder(input_dir)
    _, basename = essentials.get_path_leaf(input_dir)
    working_folder = output_dir.joinpath(basename)

    try: os.makedirs(working_folder)
    except FileExistsError: pass
    try: os.makedirs(kmc_dummy_folder)
    except FileExistsError: pass

    if sys.platform == "darwin":
        out = subprocess.run(["ulimit", "-n", "2048"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if out.returncode != 0: raise RuntimeError("Unable to increment temporary file limit on macOS")

    kmced_files = list()
    for k in ks:
        for fpath in fastxs:
            kmced_files.append(kmc_count(k, fpath, working_folder, kmc_dummy_folder, mmemory, unsorted, canonical, force))
    return kmced_files

def kmc_jaccard(executable, kmc_file1, kmc_file2) -> tuple[float, float, float]:
    executable = pathlib.Path(executable)
    kmc_file1 = pathlib.Path(kmc_file1)
    kmc_file2 = pathlib.Path(kmc_file2)
    assert executable.exists()
    # assert kmc_file1.exists()
    # assert kmc_file2.exists()

    cws_command = [executable, "full", str(kmc_file1), str(kmc_file2)]
    out = subprocess.run(cws_command, stderr=sys.stderr, stdout=subprocess.PIPE)
    if out.returncode != 0:
        raise RuntimeError("Error while performing kmc jaccard comparison\n")
    return tuple(map(float, out.stdout.decode("utf-8").split(',')))