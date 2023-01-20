#!/bin/sh

THISPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )";
REPOPATH="$(dirname -- "$THISPATH")";
KMPEXE="$REPOPATH/debug_build/kmp";

INDIR=$1;
DUMPDIR=$2; #"../../../tmp/binary_dumps"
TEMPDIR=$3; #"../../../tmp/tmp"

# echo $THISPATH;
# echo $REPOPATH;
# echo $KMPEXE;
# echo $INDIR;
# echo $DUMPDIR;
# echo $TEMPDIR;

N=10
K=15
Z=4

if [ "x$INDIR" = "x" ] ; then 
    echo "[Error] unspecified input directory where test files are located";
    exit 1 ; 
fi
if [ "x$DUMPDIR" = "x" ] ; then 
    echo "[Error] Unspecified output sketch directory";
    exit 1 ; 
fi
if [ "x$TEMPDIR" = "x" ] ; then 
    echo "[Error] Unspecified temporary directory";
    exit 1 ; 
fi

IBLT1="$DUMPDIR/BS000689.iblt.bin"
IBLT2="$DUMPDIR/BS000694.iblt.bin"

rm $IBLT1
rm $IBLT2

$KMPEXE build "$INDIR/BS000689.1.fasta" $IBLT1 -n $N -k $K -z $Z --tmp-dir $TEMPDIR
$KMPEXE build "$INDIR/BS000694.1.fasta" $IBLT2 -n $N -k $K -z $Z --tmp-dir $TEMPDIR
$KMPEXE diff $IBLT1 $IBLT2 -l "."
