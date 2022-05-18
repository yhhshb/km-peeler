function debian_realpath() {
    f=$@;
    if [ -d "$f" ]; then
        base="";
        dir="$f";
    else
        base="/$(basename "$f")";
        dir=$(dirname "$f");
    fi
    dir=$(cd "$dir" && /bin/pwd);
    echo "$dir$base";
}

function get_mutrate_from_jaccard() {
    local J=$1;
    local klen=$2;
    echo $(python3 shlib.py mfromj -j $J -k $klen);
}

function get_basename() {
    filename=$1;
    echo $filename | rev | cut -f 2- -d '.' | rev;
}

function ex5_wrapper() {
    local output=$1;
    local slen=$2;
    local klen=$3;
    local zlen=$4;
    local ntrials=$5;
    local npoints=$6;
    local inc=$(echo 4k 1 $npoints /p | dc);
    local seed=0;
    for J in $(LC_ALL=en_US.UTF-8 seq 0 $inc 1); 
    do 
        python3 $ESCRIPT ex5 -o "$output" -l $slen -k $klen -z $zlen -m $(get_mutrate_from_jaccard $J $klen) -t $ntrials --seed $seed;
        seed=$((seed + 1)); 
    done
    python3 shlib.py ex5post -i "$output" -o "$(get_basename $output).pdf";
    # for M in $(LC_ALL=en_US.UTF-8 seq 0.05 0.05 0.5); do python3 $ESCRIPT ex5 -o "$output" -l $slen -k $klen -z $zlen -m $M -t $ntrials; done
}

function ex5_alt_wrapper() {
    local output=$1;
    local slen=$2;
    local klen=$3;
    local zlen=$4;
    local ntrials=$5;
    local npoints=$6;
    local inc=$(echo 4k 1 $npoints /p | dc);
    local seed=0;
    for J in $(LC_ALL=en_US.UTF-8 seq 0 $inc 1); 
    do 
        python3 $ESCRIPT ex5alt -o "$output" -l $slen -k $klen -z $zlen -m $(get_mutrate_from_jaccard $J $klen) -t $ntrials --seed $seed;
        seed=$((seed + 1)); 
    done
}

function ex5_nojaccard_wrapper() {
    local output=$1;
    local slen=$2;
    local klen=$3;
    local zlen=$4;
    local ntrials=$5;
    local npoints=$6;
    local inc=$(echo 4k 1 $npoints /p | dc);
    local seed=0;
    for MUTR in $(LC_ALL=en_US.UTF-8 seq 0 0.01 0.26); 
    do 
        python3 $ESCRIPT ex5alt -o "$output" -l $slen -k $klen -z $zlen -m $MUTR -t $ntrials --seed $seed;
        seed=$((seed + 1)); 
    done
}

function make_out_name() {
    local prefix=$1
    local slen=$2
    local klen=$3
    local zlen=$4
    local ntrials=$5
    local ext=$6
    slen=$(python3 shlib.py intfmt -v $slen)
    klen=$(python3 shlib.py intfmt -v $klen)
    zlen=$(python3 shlib.py intfmt -v $zlen)
    ntrials=$(python3 shlib.py intfmt -v $ntrials)
    echo "${prefix}_l${slen}_k${klen}_z${zlen}_trials${ntrials}.${ext}"
}

script_name=$0
script_full_path=$(dirname $(debian_realpath "$0"))
aldiff_path=$(dirname "$script_full_path")
ESCRIPT="$script_full_path/experiments.py"
EX1POST="$script_full_path/ex1_post.py" #FIXME
RESULTS="$aldiff_path/results"
TMPFOLDER="$aldiff_path/tmp"
DFOLDER="../../../datasets"

#### EX2 #######################################################################################################
## Section Results.2 of the paper: performance measurements (syncmers only)
REX2="$RESULTS/ex2"

m=0.01
# python3 $ESCRIPT ex2 -o "$REX2/m01/l1K_k15_z4_n100_m01_e170.csv" -l 1000 -k 15 -z 4 -m $m -n 100 -e 1.7 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m01/l10K_k15_z4_n100_m01_e08.csv" -l 10000 -k 15 -z 4 -m $m -n 100 -e 0.8 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m01/l100K_k15_z4_n100_m01_e07.csv" -l 100000 -k 15 -z 4 -m $m -n 100 -e 0.7 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m01/l1M_k15_z4_n100_m01_e130.csv" -l 1000000 -k 15 -z 4 -m $m -n 100 -e 1.3 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m01/l10M_k15_z4_n10_m01_e06.csv" -l 10000000 -k 15 -z 4 -m $m -n 10 -e 0.6 --wfolder $TMPFOLDER

# python3 $ESCRIPT ex2 -o "$REX2/l100M_k15_z4_n10_m01.csv" -l 100000000 -k 15 -z 4 -m $m -n 10 --wfolder $TMPFOLDER

m=0.001
EPSILON=0.05
# python3 $ESCRIPT ex2 -o "$REX2/m001/l1K_k15_z4_n100_m001_eboh.csv" -l 1000 -k 15 -z 4 -m $m -n 100 -e 1 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m001/l10K_k15_z4_n100_m001_e60.csv" -l 10000 -k 15 -z 4 -m $m -n 100 -e 0.6 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m001/l100K_k15_z4_n100_m001_e70.csv" -l 100000 -k 15 -z 4 -m $m -n 100 -e 0.7 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m001/l1M_k15_z4_n100_m001_e75.csv" -l 1000000 -k 15 -z 4 -m $m -n 100 -e 0.75 --wfolder $TMPFOLDER
# python3 $ESCRIPT ex2 -o "$REX2/m001/l10M_k15_z4_n10_m001_e05.csv" -l 10000000 -k 15 -z 4 -m $m -n 10 -e 0.05 --wfolder $TMPFOLDER

#### EX3 ########################################################################################################
## Section Results.3: minHash vs syncmers + IBLT
REX3="$RESULTS/ex3"

COVID="$REX3/covid"
# for MEM in $(seq -f %1.0f 2000 1000 10000); do
#     python3 $ESCRIPT ex3 -i "$DFOLDER/covid50" -o "$COVID/covid50_k15_z4_q10_m$MEM.csv" -f 1 -k 15 -z 4 -c -m $MEM -w $TMPFOLDER -r 1 2 4 8 16 ;
# done

SPNEU="$REX3/spneumoniae_sican"
# for MEM in $(seq -f %1.0f 160000 20000 300000); do
#     python3 $ESCRIPT ex3 -i "$DFOLDER/spneu" -o "$SPNEU/spneumoniae_k15_z4_q10_m$MEM.csv" -f 1 -k 15 -z 4 -c -m $MEM -w $TMPFOLDER -r 1 2 4 8 16 ;
# done

#### EX5 ALTERNATE VERSION ####################################################################################
## Section Results.1 of the paper: syncmers are better than sampling and minimizers | outputs all indivdual trials for post-processing
EX5PREFIX="$RESULTS/ex5/alt"

T=500
L=10000
Z=4
# for K in 15 ; do
#     rm -f -- $(make_out_name $EX5PREFIX $L $K $Z $T "csv");
#     ex5_alt_wrapper $(make_out_name $EX5PREFIX $L $K $Z $T "csv") $L $K $Z $T 20 || exit 1;
# done

#### EX5 THIRD VERSION ####################################################################################
## 
EX5PREFIX="$RESULTS/ex5/noj"

T=500
L=10000
Z=4
for K in 15 ; do
    rm -f -- $(make_out_name $EX5PREFIX $L $K $Z $T "csv");
    ex5_nojaccard_wrapper $(make_out_name $EX5PREFIX $L $K $Z $T "csv") $L $K $Z $T 20 || exit 1;
done

#### EX6 #######################################################################################################

# EX6PREFIX="$RESULTS/ex6"
# L=1000
# Z=4
# K=15
# for M in $(LC_ALL=en_US.UTF-8 seq 0 0.01 1) ; 
# do
#     python3 $ESCRIPT ex6 -l $L -k $K -z $Z -m $M -s 42
# done

#### EX7 ######################################################################################################
REX7="$RESULTS/ex7"
# rm "$EX7PREFIX/covid_sampled.csv"
# python3 $ESCRIPT ex7 -i "$DFOLDER/covid50" -o "$REX7/covid_sampled.csv" -f 1 -k 15 -z 4 -e 0.6 -c -w $TMPFOLDER -r 1
# python3 $ESCRIPT ex7 -i "$DFOLDER/covid50" -o "$REX7/covid_sampled.csv" -f 1 -k 15 -z 4 -e 0.5 -c -w $TMPFOLDER -r 2
# python3 $ESCRIPT ex7 -i "$DFOLDER/covid50" -o "$REX7/covid_sampled.csv" -f 1 -k 15 -z 4 -e 0.5 -c -w $TMPFOLDER -r 4
# python3 $ESCRIPT ex7 -i "$DFOLDER/covid50" -o "$REX7/covid_sampled.csv" -f 1 -k 15 -z 4 -e 1.4 -c -w $TMPFOLDER -r 8
# python3 $ESCRIPT ex7 -i "$DFOLDER/covid50" -o "$REX7/covid_sampled.csv" -f 1 -k 15 -z 4 -e 2.5 -c -w $TMPFOLDER -r 16

#### EX8 ######################################################################################################
REX8="$RESULTS/ex8"
COVID="$REX8/covid"
RNDRESFOLDER="$REX8/random"
MEM=10000
# for MEM in $(seq -f %1.0f 8000 1 8000); do
# python3 $ESCRIPT ex8 -i "$DFOLDER/covid50" -o "$COVID/covid50_k15_z4_m$MEM""_c.csv" -k 15 -z 4 -c -m $MEM -w $TMPFOLDER ;
# python3 $ESCRIPT ex8 -i "$DFOLDER/covid2" -o "$COVID/covid2_k15_z4_m$MEM""_c.csv" -k 15 -z 4 -m $MEM -w $TMPFOLDER ;
# done


RND50M01="$DFOLDER/random_n50_l30k_m01_S0"
RND50M01OUT="$RNDRESFOLDER/random_n50_l30k_m01_S0_k15_z4_c_size$MEM.csv"
RND50M001="$DFOLDER/random_n50_l30k_m001_S1"
RND50M001OUT="$RNDRESFOLDER/random_n50_l30k_m001_S0_k15_z4_c_size$MEM.csv"
# mkdir $RND50M01
MEM=150000
# python3 $ESCRIPT rndset -l 30000 -n 49 -m 0.01 -s 0 -o $RND50M01
# python3 $ESCRIPT ex8 -i $RND50M01 -o $RND50M01OUT -k 15 -z 4 -c -m $MEM -w $TMPFOLDER ;
# mkdir $RND50M001
# MEM=100000
# python3 $ESCRIPT rndset -l 30000 -n 49 -m 0.001 -s 1 -o $RND50M001
# python3 $ESCRIPT ex8 -i $RND50M001 -o $RND50M001OUT -k 15 -z 4 -c -m $MEM -w $TMPFOLDER ;

## -----------------------------------------------------------------------------------------------------------------------------------------

# python3 $ESCRIPT ex8 -i "$DFOLDER/covid50" -o "$COVID/covid50_k15_z4_m$MEM""_x.csv" -k 15 -z 4 -x -m $MEM -w $TMPFOLDER ;
# python3 $ESCRIPT ex8 -i "$DFOLDER/covid50" -o "$COVID/covid50_k15_z4_m$MEM.csv" -k 15 -z 4 -m $MEM -w $TMPFOLDER ;
# python3 $ESCRIPT ex8 -i "$DFOLDER/covid50" -o "$COVID/covid50_k15_z4_m$MEM""_c_x.csv" -k 15 -z 4 -c -x -m $MEM -w $TMPFOLDER ;

#### EX1 #####################################################################################################
## NOT USED in the paper: 
## finding true "n" parameter for IBF construction on random sequences
REX1="$RESULTS/ex1"

m=0.01
# python3 $ESCRIPT ex1 -o "$REX1/l1k_k15_z4_n1000.csv" -l 1000 -k 15 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l1k_k15_z8_n1000.csv" -l 1000 -k 15 -z 8 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l1k_k12_z4_n1000.csv" -l 1000 -k 12 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l1k_k10_z4_n1000.csv" -l 1000 -k 10 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l1k_k8_z4_n1000.csv" -l 1000 -k 8 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l10k_k15_z4_n1000.csv" -l 10000 -k 15 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l10k_k15_z8_n1000.csv" -l 10000 -k 15 -z 8 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l10k_k12_z4_n1000.csv" -l 10000 -k 12 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l10k_k10_z4_n1000.csv" -l 10000 -k 10 -z 4 -m $m -n 1000
# python3 $ESCRIPT ex1 -o "$REX1/l10k_k8_z4_n1000.csv" -l 10000 -k 8 -z 4 -m $m -n 1000

# python3 $EX1POST "$REX1/l1k_k15_z4_n1000.csv" "$REX1/l1k_k15_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l1k_k15_z8_n1000.csv" "$REX1/l1k_k15_z8_n1000.pdf"
# python3 $EX1POST "$REX1/l1k_k12_z4_n1000.csv" "$REX1/l1k_k12_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l1k_k10_z4_n1000.csv" "$REX1/l1k_k10_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l1k_k8_z4_n1000.csv" "$REX1/l1k_k8_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l10k_k15_z4_n1000.csv" "$REX1/l10k_k15_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l10k_k15_z8_n1000.csv" "$REX1/l10k_k15_z8_n1000.pdf"
# python3 $EX1POST "$REX1/l10k_k12_z4_n1000.csv" "$REX1/l10k_k12_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l10k_k10_z4_n1000.csv" "$REX1/l10k_k10_z4_n1000.pdf"
# python3 $EX1POST "$REX1/l10k_k8_z4_n1000.csv" "$REX1/l10k_k8_z4_n1000.pdf"

#### EX5 ########################################################################################################
## NOT USED in the paper: 
## Section Results.1 of the paper: syncmers are better than sampling and minimizers | include variance in the csv
EX5PREFIX="$RESULTS/ex5/runs"

T=50
L=10000

# Shown results
Z=4
# for K in 15 127; do
#     rm -f -- $(make_out_name $EX5PREFIX $L $K $Z $T "csv");
#     produce_results $(make_out_name $EX5PREFIX $L $K $Z $T "csv") $L $K $Z $T 50 || exit 1;
# done

# Control
# Z=13
# for K in 15 ; do
#     rm -f -- $(make_out_name $EX5PREFIX $L $K $Z $T "csv");
#     produce_results $(make_out_name $EX5PREFIX $L $K $Z $T "csv") $L $K $Z $T 50 || exit 1;
# done

SPNEU="$REX3/spneumoniae_nocan"
## NOT GOOD!!!
# for MEM in $(seq -f %1.0f 1500000 100000 1700000); do
#     python3 $ESCRIPT ex3 -i "$DFOLDER/spneu" -o "$SPNEU/spneumoniae_k15_z4_q10_m$MEM.csv" -f 1 -k 15 -z 4 -m $MEM -w $TMPFOLDER ;
# done