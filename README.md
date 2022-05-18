## <a name="uguide"></a>Users' Guide

km-peeler is an implementation of Invertible Bloom Lookup Tables [1][belbasi] targeting k-mer sets.
IBLTs are a generalization of Bloom filters for storing sets of elements (keys).
In addition to insertions, IBLTs also support deletions and listing of the keys they contain.
The latter operation is based on the peelability properties of random hypergraphs [2][pagh] and it succeeds, with high probability, if the number of items in an IBLT is below a given threshold used during bulding.
This condition must be satisfied only when items have to be retrieved and a overly-charged IBLT can return to be peelable if enough items are deleted.

The above property can be applied to quickly find symmetric differences of k-mer sets in space that only depends on the size of the difference and not on the size of the involved sets.

## <a name="install"></a>Installation

km-peeler is written in C99 and was tested on both macOS and Linux.
All dependencies for the main program are already included.

On the other hand, test/plotting scripts are written in python 3 (with some helper tools in C++) and need the following libraries to work properly:
- psutil
- pandas
- matplotlib
- seaborn

Installation is done by cloning this repository and running the configuration script *configure.py* (which calls make internally)

```sh
git clone https://github.com/yhhshb/km-peeler.git
cd km-peeler
python3 configure.py <type> <options>
```

Calling configure.py generates a header file containing the maximum length of the payload of each bucket.
Examples:

### Configuring km-peeler for fixed-length sequences (k-mers)

```sh
python3 configure.py sequences -a <algorithm> -k <sequence length (number of bases)> -z <z length (number of bases)>
```

configures km-peeler for fixed-length sequence storage.
**algorithm** can be one of the following:
- kmers --> option <z> is ignored
- syncmers --> option <z> (< <k>) is the length of the z-mers used to define syncmers
- minimizers --> option <z> is the length of a k-mer's length. This configure km-peeler for storing pairs of k-mers grouped by their minimizers (so sequences of length 2*k-z). Normal minimizers are equivalent to k-mers or syncmers since they have length =k.

### Configuring km-peeler for hashes

```sh
python3 configure.py hashes -l <hash length in bits> -s <number of minimum hashes>
```

In particular, for 64-bit hashes:
```sh
python3 configure.py hashes -l 64
```

and minHash sketches made of s=1000, 64-bit hash values:

```sh
python3 configure.py hashes -l 64 -s 1000
```
Note that minHash sketches are seen as a very big hash value.

### Configuring km-peeler for variable-length sequences

```sh
python3 configure.py vla -l <maximum sequence length (number of bases)>
```

## <a name="general"></a>General usage

The main executable of km-peeler is called ibltseq with its sub-command `build` expecting a **set** of strings in input.
Multi-sets of strings can be obtained by calling ibltseq with the subcommands [kmers, minimizers, syncmers].
A simple python script (2set.py) is provided to remove duplicates.

For example, in order to constructing an IBLT on the set of k-mers of a sequence is done by running:
```sh
ibltseq kmers -k 15 -i <input.fasta> | python3 2set.py | ibltseq build -n <IBLT threshold> -o <output IBLT>
```

[belbasi]: https://doi.org/10.48550/arXiv.1101.2245
[pagh]: https://doi.org/10.1007/978-3-642-14165-2_19
