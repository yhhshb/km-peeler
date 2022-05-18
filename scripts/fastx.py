import gzip

enc_nuc_table = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

"""
enc_nuc_table_extended deals with letters different from 
{A, C, G, T, N} by simply replacing them with one non-ambigous DNA
U -> T
R -> G | Y -> C | K -> G | M -> A | S -> C |W -> A
B -> A | D -> C | H -> G | V -> T
Note that B, D, H and V are mapped to the letters they
do not represent for sure. This breaks the definition
but it is only to have a deterministic choice for the
replacement. DO NOT USE THIS TABLE other than set comparisons.
"""
enc_nuc_table_extended = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4, #5 is the gap character '-' -> it should be filtered out
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 0, 1,  1, 4, 4, 2,  2, 4, 4, 2,  4, 0, 4, 4, 
    4, 4, 3, 2,  3, 3, 3, 0,  4, 2, 4, 4,  4, 4, 4, 4, 
    4, 0, 0, 1,  1, 4, 4, 2,  2, 4, 4, 2,  4, 0, 4, 4, 
    4, 4, 3, 2,  3, 3, 3, 4,  4, 2, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

dec_nuc_table = ['A', 'C', 'G', 'T', 'N']

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def is_gzipped(fastx_path: str) -> bool:
    gzipped = True
    with gzip.open(fastx_path, 'r') as fh:
        try:
            fh.read(1)
        except OSError:
            gzipped = False
    return gzipped

def read(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l.rstrip() # save this line #FIXED removed l[:-1] because last char was ignored
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l.rstrip()#FIXED removed l[:-1] because last char was ignored
                break
            seqs.append(l.rstrip())#FIXED removed l[:-1] because last char was ignored
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l.rstrip())#FIXED removed l[:-1] because last char was ignored
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def fasta_write(fp, name: str, seq: str, last: bool = True):
    if last: fp.write(">{}\n{}".format(name, seq))
    else: fp.write(">{}\n{}\n".format(name, seq))

def length(fp) -> int:
    """Return sum of contig lengths"""
    total_length = 0
    for _, seq, _ in read(fp):
        total_length += len(seq)
    return total_length
