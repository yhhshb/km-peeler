import os
import sys
import pathlib
import string
import random

gzipped_extensions = ["gz"]
fasta_extensions = ["fasta", "fa", "fna"]
fastq_extensions = ["fastq", "fq", "fnq"]

def get_random_name(extension: str, size=20, chars=string.ascii_uppercase + string.digits) -> str:
    suffix = ('.' + extension) if extension else ""
    return "".join(random.choice(chars) for _ in range(size)) + suffix

def make_txt_to_fasta(txt: str, fasta: str):
    counter = 0
    with open(txt, "r") as sh:
        with open(fasta, "w") as fh:
            for line in sh:
                fh.write(">{}\n{}\n".format(counter, line.strip())) #strip because line has \n at the end
                counter += 1

def is_empty_file(file: str):
    try:
        with open(file, "r") as csvh:
            csvh.seek(0, os.SEEK_END)
            if not csvh.tell(): write_header = True
            else: write_header = False
    except FileNotFoundError as e:
        write_header = True
    return write_header

def get_path_leaf(path) -> tuple[str]: # path can also be of type Path (no method strip)
    path = pathlib.Path(path)
    parent = path.parent
    node = path.name
    return parent, node

def get_fastx_format(filename):
    exts = pathlib.Path(filename).suffixes
    file_type = None
    gzipped = False
    for suf in exts:
        suf = suf.strip('.')
        if suf.lower() in fasta_extensions:
            if file_type: raise ValueError("Multiple suitable fastx extensions")
            file_type = "fasta"
        elif suf.lower() in fastq_extensions:
            if file_type: raise ValueError("Multiple suitable fastx extensions")
            file_type = "fastq"
        elif (suf in gzipped_extensions) and (suf == exts[-1]):
            gzipped = True
    if file_type == None: raise ValueError("Not in fastx format")
    return (file_type, gzipped)

def get_fastx_in_folder(path) -> list[str]:
    path = pathlib.Path(path)
    fastx_files = list()
    for f in path.iterdir():
        if f.is_file():
            try:
                get_fastx_format(f)
                fastx_files.append(str(f))
            except ValueError as e:
                sys.stderr.write("Failed to categorise file {}: {}\n".format(f, e))
    fastx_files.sort(key = lambda v : v.upper())
    return fastx_files