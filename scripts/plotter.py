import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def jaccard_postprocessing_main(args):
    dfs = [pd.DataFrame(pd.read_csv(file)) for file in args.i]
    all = pd.concat(dfs)
    dfs = None
    xname = "Allocated space (Bytes)"
    yname = "Average absolute error"
    # all.dropna(inplace=True)
    #all["mhj"] = all["mhj"].round(decimals=3)
    all = all[all["query"] != all["reference"]]
    all["mherr"] = (all["mhj"] - all["exj"]).abs()
    # all["syncmherr"] = (all["syncmhj"] - all["exj"]).abs()

    sampling_rates = list()
    for col in all.columns:
        if "syncj r=" in col:
            rate = int(col.replace("syncj r=", ""))
            sampling_rates.append(rate)

    syncmers_columns = list()
    for rate in sampling_rates:
        abs_err_col = r"syncmers $\eta$={}".format(rate)
        iblt_col = "syncj r={}".format(rate)
        all[iblt_col] = pd.to_numeric(all[iblt_col], errors="coerce")
        all[abs_err_col] = (all[iblt_col] - all["exj"]).abs()
        syncmers_columns.append(abs_err_col)
    nandf = all[syncmers_columns].rename(columns={"size":xname}).isnull()
    nandf[xname] = all["size"]
    all.dropna(inplace=True)
    # errdf = all[["size", "mherr", "syncmherr"] + syncmers_columns]
    errdf = all[["size", "mherr"] + syncmers_columns]
    df = errdf.groupby(["size"], as_index=False).mean()
    # df.rename(columns={"size":xname, "mherr":"minHash", "syncmherr":"syncmers + minHash"}, inplace=True)
    df.rename(columns={"size":xname, "mherr":"minHash"}, inplace=True)
    melted_err = df.melt(xname, var_name="method", value_name=yname)
    fig, ax = plt.subplots()
    sns.barplot(data=melted_err, x=xname, y=yname, hue="method", ax=ax)
    lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) # legend outside plot area
    fig.savefig(args.o + "_errors.pdf", bbox_inches="tight", bbox_extra_artists=(lgd,))
    nandf = nandf.groupby([xname], as_index=False).sum().astype(int)#.reset_index(name="count")
    melted_nandf = nandf.melt(xname, var_name="method", value_name="Nan count")
    fig, ax = plt.subplots()
    sns.barplot(data=melted_nandf, x=xname, y="Nan count", hue="method", ax=ax)
    fig.savefig(args.o + "_nan.pdf", bbox_inches="tight")

def main(args):
    if args.command == "jaccard": return jaccard_postprocessing_main(args)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_jaccard = subparsers.add_parser("jaccard", help="Post-processing for kmp jaccard experiments")
    parser_jaccard.add_argument("-i", help="list of csv files produced by ex3 for different pre-allocated space", type=str, nargs='+', required=True)
    parser_jaccard.add_argument("-o", help="name prefix for the plots", type=str, required=True)

    return parser

if __name__ == "__main__":
    parser = parser_init()
    args = parser.parse_args(sys.argv)
    main(args)