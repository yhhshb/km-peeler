import os
import sys
import math
import shutil
from numpy import require
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

format_letters = ['', 'K', 'M', 'G', 'T']
format_time = ['n', r'$\micro$', 'm', '']
fmt_to_mag = {'':0,'K':1,'M':2,'G':3,'T':4}

def get_time_magnitude(ns: int):
    magnitude = 0
    while(magnitude < 3 and ns > 1000):
        magnitude += 1
        ns /= 1000
    return magnitude

def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), format_letters[magnitude])

def intfmt_main(args):
    sys.stdout.write(human_format(args.v))

def mfromj_main(args):
    try:
        mut_rate = -1/args.k * math.log(2*args.j/(1+args.j))
    except ValueError:
        mut_rate = 1
    sys.stdout.write("{}".format(mut_rate))

def ffilter_main(args):
    assert os.path.isfile(args.l)
    assert os.path.isdir(args.i)
    assert os.path.isdir(args.o)
    names = set()
    with open(args.l, "r") as lh:
        for name in lh:
            names.add(os.path.basename(os.path.splitext(name.strip())[0]))
    file_ls = [f for f in os.listdir(args.i) if os.path.isfile(os.path.join(args.i, f))]
    file_ls = [f for f in file_ls if os.path.splitext(f)[0] in names]
    for f in file_ls:
        full_ipath = os.path.join(args.i, f)
        full_opath = os.path.join(args.o, f)
        shutil.copyfile(full_ipath, full_opath)

def ex1post_main(args):
    df = pd.DataFrame(pd.read_csv(args.i))
    df["diff"] = df["nmm"] - df["nsm"]
    mxdiff = max(df["diff"])
    fig, ax = plt.subplots()
    ax = df.hist(column="diff", grid=False, bins=mxdiff, ax=ax)
    plt.xlabel("minimum n for minimizers - minimum n for syncmers")
    plt.ylabel("# of values")
    plt.title("Minimum n_mm - n_sm for which the sketch becomes peelable")
    fig.savefig(args.o)

def ex2post_main(args):
    dfs = [pd.DataFrame(pd.read_csv(file)) for file in args.i]
    all = pd.concat(dfs)
    dfs = None # possible deallocation
    # plot space comparison
    df = all[all["peelable"] == True].groupby(["length", "k", "z", "mutp", "indelf", "extp"], as_index=False).mean()
    df["kmcsize"] = df["kmcosize"] + df["kmcmsize"] #size of original and mutated version, column names are not significant
    df["ibfsize"] = df["ibfsize"] * 2
    diskdf = df[["length", "ibfsize", "kmcsize"]].copy()
    if args.smagnitude: magnitude = fmt_to_mag[args.smagnitude]
    else: magnitude = int(diskdf[["ibfsize", "kmcsize"]].min().min() / 1000)
    print(diskdf[["ibfsize", "kmcsize"]].min().min())
    print(magnitude)
    diskdf["ibfsize"] = diskdf["ibfsize"] / (1000**magnitude)
    diskdf["kmcsize"] = diskdf["kmcsize"] / (1000**magnitude)
    diskdf.rename(columns={"length":"L", "ibfsize":"IBLT", "kmcsize":"kmc"}, inplace=True)
    yname = "Disk space ({}B)".format(format_letters[magnitude])
    meltedisk = diskdf.melt("L", var_name="method", value_name=yname)
    fig, ax = plt.subplots()
    sns.barplot(data=meltedisk, x="L", y=yname, hue="method", ax=ax)
    fig.savefig(args.o + "_space.pdf", bbox_inches="tight")
    # plot time comparisons
    df["ibftime"] = df["obtime"] + df["mbtime"]
    df["kmctime"] = df["kmcotime"] + df["kmcmtime"]
    timedf = df[["length", "ibftime", "kmctime"]].copy()
    #print(timedf[["ibftime", "kmctime"]].min().min())
    magnitude = get_time_magnitude(timedf[["ibftime", "kmctime"]].min().min())
    #print(magnitude)
    timedf["ibftime"] = timedf["ibftime"] / (1000**magnitude)
    timedf["kmctime"] = timedf["kmctime"] / (1000**magnitude)
    timedf.rename(columns={"length":"L", "ibftime":"IBLT", "kmctime":"kmc"}, inplace=True)
    yname = "construction time ({}s)".format(format_time[magnitude])
    meltedtime = timedf.melt("L", var_name="method", value_name=yname)
    fig, ax = plt.subplots()
    sns.barplot(data=meltedtime, x="L", y=yname, hue="method", ax=ax)
    fig.savefig(args.o + "_ctime.pdf", bbox_inches="tight")

def ex3post_main(args):
    dfs = [pd.DataFrame(pd.read_csv(file)) for file in args.i]
    all = pd.concat(dfs)
    dfs = None
    xname = "Allocated space (B)"
    yname = "Average absolute error"
    # all.dropna(inplace=True)
    all["mhj"] = all["mhj"].round(decimals=3)
    all["mherr"] = (all["mhj"] - all["exj"]).abs()
    all["syncmherr"] = (all["syncmhj"] - all["exj"]).abs()
    sampling_rates = list()
    for col in all.columns:
        if "syncj r=" in col:
            rate = int(col.replace("syncj r=", ""))
            sampling_rates.append(rate)
    syncmers_columns = list()
    for rate in sampling_rates:
        col = r"syncmers $\nu$={}".format(rate)
        all[col] = (all["syncj r={}".format(rate)] - all["exj"]).abs()
        syncmers_columns.append(col)
    # all["syncmers"] = (all["syncj"] - all["exj"]).abs()
    # all["smplsyncerr"] = (all["smplsyncj"] - all["exj"]).abs()
    # rename(columns={"size":xname})
    nandf = all[syncmers_columns].rename(columns={"size":xname}).isnull()
    nandf[xname] = all["size"]
    all.dropna(inplace=True)
    errdf = all[["size", "mherr", "syncmherr"] + syncmers_columns]
    df = errdf.groupby(["size"], as_index=False).mean()
    df.rename(columns={"size":xname, "mherr":"minHash", "syncmherr":"syncmers + minHash"}, inplace=True)
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

    # df = pd.DataFrame(pd.read_csv(args.i))
    # print(df.corr("spearman"))# .style.background_gradient(cmap="Blues")
    # mh_sum_of_errors = (abs(df["exj"] - df["mhj"])).sum() / df.shape[0]
    # # smpl_sum_of_errors = (abs(df["exj"] - df["smplj"])).sum() / df.shape[0]
    # sync_sum_of_errors = (abs(df["exj"] - df["syncj"])).sum() / df.shape[0]
    # mh_sum_of_sq_errors = ((df["exj"] - df["mhj"])**2).sum() / df.shape[0]
    # # smpl_sum_of_sq_errors = ((df["exj"] - df["smplj"])**2).sum() / df.shape[0]
    # sync_sum_of_sq_errors = ((df["exj"] - df["syncj"])**2).sum() / df.shape[0]
    # print("minHash mean error: absolute = {0:.6f} squared = {0:.6f}".format(mh_sum_of_errors, mh_sum_of_sq_errors))
    # # print("sampling mean error: absolute = {0:.6f} squared = {0:.6f}".format(smpl_sum_of_errors, smpl_sum_of_sq_errors))
    # print("syncmers mean error: absolute = {0:.6f} squared = {0:.6f}".format(sync_sum_of_errors, sync_sum_of_sq_errors))

def ex5post_main(args):
    df = pd.DataFrame(pd.read_csv(args.i))
    fig, ax = plt.subplots()
    clrs = sns.color_palette()
    ax.plot([0, 1], [0, 1], linewidth=1, color=clrs[3])
    ax.scatter(df["exact jaccard"], df["sampled jaccard"], marker='.', s=8, color="black", label="sampling")
    ax.scatter(df["exact jaccard"], df["syncmers jaccard"], marker='+', s=8, color="black", label="syncmers")
    ax.scatter(df["exact jaccard"], df["minimizers jaccard"], marker='x', s=8, color="black", label="minimizers")
    ax.set_xlabel(r"$\mathtt{J}$") # \: for space
    ax.set_ylabel(r"$\widehat{\mathtt{J}}$", rotation=0, ha="right")
    ax.legend()
    # ax = sns.scatterplot(data=df, x="exact jaccard", y="syncmers jaccard", s=10, markers=["x", "o"])
    # sns.scatterplot(data=df, x="exact jaccard", y="minimizers jaccard", s=10, markers=['x'], ax=ax)
    # sns.lineplot(data=df, x="exact jaccard", y="exact jaccard", color='r', linewidth=0.1, ax=ax)
    # fig = ax.get_figure()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    fig.savefig(args.o, bbox_inches="tight")
    
    df["var diff"] = df["var sampled"] - df["var syncmers"]
    fig, ax = plt.subplots()
    g2 = sns.barplot(x="exact jaccard", y="var diff", color="black", data=df, ax=ax)
    # ax.scatter(df["exact jaccard"], df["var diff"], marker='v', color="black")
    g2.set(xticklabels=[])  
    # g2.set(title='Penguins: Body Mass by Species for Gender')
    g2.set(xlabel=None)
    g2.tick_params(bottom=False)  # remove the ticks
    ax.set_xlabel(None) # \: for space
    ax.set_ylabel(r"$\sigma_{ns}^{2} - \sigma_{s}^{2}$", ha="right")
    # samplevar = df[["exact jaccard", "var sampled"]].rename(columns={"var sampled":"variance"})
    # samplevar["method"] = "sampling"
    # syncvar = df[["exact jaccard", "var syncmers"]].rename(columns={"var syncmers":"variance"})
    # syncvar["method"] = "syncmers"
    # variances = pd.concat([samplevar, syncvar])
    # fg = sns.catplot(x="exact jaccard", y="variance", hue="method", data=variances, kind="bar", legend_out=False)
    basename = os.path.splitext(args.o)[0]
    fig.savefig(basename + "_variance.pdf", bbox_inches="tight")

def ex5altpost_main(args):
    df = pd.DataFrame(pd.read_csv(args.i))
    clrs = sns.color_palette()
    # ------------- print same plot as ex5post
    av = df[["mutp", "exact jaccard", "sampled jaccard", "minimizers jaccard", "syncmers jaccard"]]
    av = av.groupby(["mutp"], as_index=False).mean()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], linewidth=1, color=clrs[3])
    ax.scatter(av["exact jaccard"], av["sampled jaccard"], marker='.', s=8, color="black", label="sampling")
    ax.scatter(av["exact jaccard"], av["syncmers jaccard"], marker='+', s=8, color="black", label="syncmers")
    ax.scatter(av["exact jaccard"], av["minimizers jaccard"], marker='x', s=8, color="black", label="minimizers")
    ax.set_xlabel(r"$\mathtt{J}$") # \: for space
    ax.set_ylabel(r"$\widehat{\mathtt{J}}$", rotation=0, ha="right")
    ax.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    fig.savefig(args.o + "_bias.pdf", bbox_inches="tight")
    av = None
    # ------------- print sampling and syncmers errors in the same plot
    df["sampling error"] = (df["sampled jaccard"] - df["exact jaccard"]).abs()
    df["syncmers error"] = (df["syncmers jaccard"] - df["exact jaccard"]).abs()
    av = df.groupby(["mutp"], as_index=False).mean()
    fig, ax = plt.subplots()
    ax.plot(av["exact jaccard"], av["syncmers error"], marker='+', markersize=8, color=clrs[1], label="syncmers")
    ax.plot(av["exact jaccard"], av["sampling error"], marker='.', markersize=8, color=clrs[0], label="sampling")
    ax.set_xlabel(r"$\mathtt{J}$")
    ax.set_ylabel("Average absolute error")
    ax.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    fig.savefig(args.o + "_errors.pdf", bbox_inches="tight")
    # ------------- same thing as before but with mutation rate as the x-axis
    fig, ax = plt.subplots()
    ax.plot(av["mutp"], av["syncmers error"], marker='+', markersize=8, color=clrs[1], label="syncmers")
    ax.plot(av["mutp"], av["sampling error"], marker='.', markersize=8, color=clrs[0], label="sampling")
    ax.set_xlabel("Mutation probability")
    ax.set_ylabel("Average absolute error")
    ax.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    fig.savefig(args.o + "_mut_errors.pdf", bbox_inches="tight")
    # ------------- print point cloud
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], linewidth=1, color=clrs[3])
    ax.scatter(df["exact jaccard"], df["syncmers jaccard"], marker='o', s=0.1, color=clrs[1], label="syncmers")
    ax.scatter(df["exact jaccard"], df["sampled jaccard"], marker='o', s=0.1, color=clrs[0], label="sampling")
    ax.set_xlabel(r"$\mathtt{J}$") # \: for space
    ax.set_ylabel(r"$\overline{\mathtt{J}}$", rotation=0, ha="right")
    ax.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    fig.savefig(args.o + "_cloud.pdf", bbox_inches="tight")

def ex7post_main(args):
    df = pd.DataFrame(pd.read_csv(args.i))
    xname = r"$\nu$"
    yname = "Disk space (B)"
    df.dropna(inplace=True)
    datadf = df[["sampling rate", "size", "exj", "syncj"]]
    gb = datadf.groupby(["sampling rate"], as_index=False).mean()
    gb.rename(columns={"sampling rate":xname, "size":yname}, inplace=True)
    # melted_err = df.melt(xname, var_name="method", value_name=yname)
    fig, ax = plt.subplots()
    sns.barplot(data=gb, x=xname, y=yname, color="black", ax=ax)
    #plt.show()
    fig.savefig(args.o + "_sampled_size.pdf", bbox_inches="tight")

# def ex8post_main(args):
#     ex3_df = pd.DataFrame(pd.read_csv(args.ex3output))
#     ex8_df = pd.DataFrame(pd.read_csv(args.ex8output))
#     exj_df = ex3_df[["reference", "query", "exj"]]
#     ex8_df = pd.merge(exj_df, ex8_df, on=["reference", "query"], how="inner")
#     err3 = (ex3_df["syncj r=1"] - ex3_df["exj"]).abs().mean()
#     err8 = (ex8_df["esyncj"] - ex8_df["exj"]).abs().mean()
#     _, ax = plt.subplots()
#     sns.barplot(x=["full syncmers", "extended syncmers"], y=[err3, err8], ax=ax)
#     plt.show()

def ex8post_main(args):
    ex8_df = pd.DataFrame(pd.read_csv(args.ex8output))
    ex8_df["difference"] = (ex8_df["estimated symmetric size"] - ex8_df["true symmetric size"])
    # err8 = err8_df.mean()
    print(ex8_df[["true symmetric size", "estimated symmetric size", "difference"]].describe())
    sns.displot(data=ex8_df, x="difference", binwidth=1)
    plt.show()

def main(args):
    if (args.command == "intfmt"): return intfmt_main(args)
    elif (args.command == "mfromj"): return mfromj_main(args)
    elif (args.command == "ffilter"): return ffilter_main(args)
    elif (args.command == "ex1post"): return ex1post_main(args)
    elif (args.command == "ex2post"): return ex2post_main(args)
    elif (args.command == "ex3post"): return ex3post_main(args)
    elif (args.command == "ex5post"): return ex5post_main(args)
    elif (args.command == "ex5altpost"): return ex5altpost_main(args)
    elif (args.command == "ex7post"): return ex7post_main(args)
    elif (args.command == "ex8post"): return ex8post_main(args)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_intfmt = subparsers.add_parser("intfmt", help="Format integer numbers for humans, i.e. add K, M, B, T")
    parser_intfmt.add_argument("-v", help="integer value", type=int, required=True)

    parser_mfromj = subparsers.add_parser("mfromj", help="Get mutation rate to obtain given jaccard")
    parser_mfromj.add_argument("-j", help="Jaccard value", type=float, required=True)
    parser_mfromj.add_argument("-k", help="k-mer length", type=int, required=True)

    parser_ffilter = subparsers.add_parser("ffilter", help="copy files from one folder to another based on list of names")
    parser_ffilter.add_argument("-i", help="input folder", type=str, required=True)
    parser_ffilter.add_argument("-o", help="output folder", type=str, required=True)
    parser_ffilter.add_argument("-l", help="list of names", type=str, required=True)

    parser_ex1post = subparsers.add_parser("ex1post", help="Post-processing for ex1")
    parser_ex1post.add_argument("-i", help="csv file produced by ex1", type=str, required=True)
    parser_ex1post.add_argument("-o", help="plot img file", type=str, required=True)

    parser_ex2post = subparsers.add_parser("ex2post", help="Post-processing for ex2")
    parser_ex2post.add_argument("-i", help="csv files produced by ex2" ,type=str, nargs='+', required=True)
    parser_ex2post.add_argument("-o", help="name prefix for the plots", type=str, required=True)
    parser_ex2post.add_argument("--smagnitude", help="space magnitude used to rescale all columns", type=str, default="", choices=format_letters)

    parser_ex3post = subparsers.add_parser("ex3post", help="Post-processing for ex3")
    parser_ex3post.add_argument("-i", help="csv files produced by ex3 for different pre-allocated space", type=str, nargs='+', required=True)
    parser_ex3post.add_argument("-o", help="name prefix for the plots", type=str, required=True)

    parser_ex5post = subparsers.add_parser("ex5post", help="Post-processing for ex5(a.k.a. plot csv table)")
    parser_ex5post.add_argument("-i", help="csv file produced by ex5", type=str, required=True)
    parser_ex5post.add_argument("-o", help="plot img file", type=str, required=True)

    parser_ex5post = subparsers.add_parser("ex5altpost", help="Post-processing for ex5 alternative version (a.k.a. plot csv table)")
    parser_ex5post.add_argument("-i", help="csv file produced by ex5 alternative version", type=str, required=True)
    parser_ex5post.add_argument("-o", help="prefix for plot files. Generated plots are: TODO", type=str, required=True)

    parser_ex7post = subparsers.add_parser("ex7post", help="Post-processing for ex7")
    parser_ex7post.add_argument("-i", help="csv file produced by ex7", type=str, required=True)
    parser_ex7post.add_argument("-o", help="name prefix for the plots", type=str, required=True)

    parser_ex8post = subparsers.add_parser("ex8post", help="Post-processing for ex8")
    # parser_ex8post.add_argument("--ex3output", help="csv file produced by ex3", type=str, required=True)
    parser_ex8post.add_argument("-i", "--ex8output", help="csv file produced by ex8", type=str, required=True)

    return parser

if __name__ == "__main__":
    parser = parser_init()
    args = parser.parse_args(sys.argv)
    main(args)