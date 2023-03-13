import sys
import pathlib
import shutil

'''
Given a "mash triangle -E" output, return a cluster of sequences that are close to each other.
'''

class MashEntry:
    def __init__(self, entry_line: str):
        el = entry_line.strip().split('\t')
        self.seq1 = el[0]
        self.seq2 = el[1]
        self.dist = float(el[2])
        self.pval = float(el[3])
        self.intersection_size, self.union_size = map(int, el[4].split('/'))
        self.jaccard = self.intersection_size / self.union_size
    
    def __str__(self) -> str:
        return "{}\t{}\t{}\t{}\t{}/{}={}".format(self.seq1, self.seq2, self.dist, self.pval, self.intersection_size, self.union_size, self.jaccard)

if __name__ == "__main__":
    entries = list()
    with open(sys.argv[1], "r") as mash_triangle_h:
        for line in mash_triangle_h:
            entries.append(MashEntry(line))
    entries.sort(key=lambda e: e.dist)
    jaccard_cutoff = float(sys.argv[2])
    cluster = set()
    for entry in entries: 
        if entry.seq1 != entry.seq2 and entry.jaccard > jaccard_cutoff:
            cluster.add(entry.seq1)
            cluster.add(entry.seq2)
    sys.stderr.write("Found cluster of {} sequences\n".format(len(cluster)))
    out_dir = pathlib.Path(sys.argv[3])
    if not out_dir.exists(): sys.stderr.write("Silent mode, the output directory does not exists\n")
    for seq in cluster:
        fasta_file = pathlib.Path(seq)
        assert fasta_file.exists()
        out_file = out_dir.joinpath(fasta_file.name)
        if out_dir.exists(): shutil.copy(fasta_file, out_file)
    for entry in entries:
        if entry.seq1 in cluster and entry.seq2 in cluster:
            print(entry)
    