import os
import sys
import pathlib
import subprocess

salmonella_url = "https://zenodo.org/record/4338293/files/assemblies.zip"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python {} <output folder>\n".format(sys.argv[0]))
        raise RuntimeError("Too few arguments")
    dataset_folder = pathlib.Path(sys.argv[1])
    print("Creating new dataset folders in {}".format(dataset_folder))
    if not dataset_folder.exists():
        os.makedirs(dataset_folder)
    # salmonella_dataset_folder = dataset_folder.joinpath("salmonella")
    # if not salmonella_dataset_folder.exists():
    #     os.makedirs(salmonella_dataset_folder)
    out = subprocess.run(["wget", salmonella_url], stderr=sys.stderr, stdout=subprocess.DEVNULL)
    if out.returncode != 0:
        raise RuntimeError("Unable to retrieve the salmonella dataset")
    salmonella_zip_file = pathlib.Path("assemblies.zip")
    assert salmonella_zip_file.exists()
    salmonella_zip_file.rename(dataset_folder.joinpath(salmonella_zip_file))
    
    tuberculosis_dataset_folder = dataset_folder.joinpath("tuberculosis")
    if not tuberculosis_dataset_folder.exists():
        os.makedirs(tuberculosis_dataset_folder)
    tuberculosis_accessions = pathlib.Path("tuberculosis_accessions.txt")
    assert tuberculosis_accessions.exists()
    with open(tuberculosis_accessions, "r") as th:
        for accession in th:
            accession = accession.strip()
            subprocess.run(["ncbi-acc-download", "--format", "fasta", accession], stderr=sys.stderr, stdout=subprocess.DEVNULL)
            downloaded = pathlib.Path(accession + ".fa") # Do not use with_suffix since names can be X.1
            assert downloaded.exists()
            downloaded.rename(tuberculosis_dataset_folder.joinpath(downloaded))
