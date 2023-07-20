#!/usr/bin/env python
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import argparse


def add_patch(legend, fc, label):
    from matplotlib.patches import Patch

    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    for i in range(len(fc)):
        handles.append(Patch(facecolor=fc[i], edgecolor=fc[i]))
        labels.append(label[i])

    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())
    return legend


def lineplots(file):
    t = Table.read(file, format="ipac")
    colors = ["red", "tan", "lime"]
    df = t.to_pandas()  # .sort_values('N_rows')
    time_cols = ["Load", "Index", "Create_Obsblocks", "Baseline_1", "Baseline_2", "Baseline_3"]
    tc = [s.replace("_", " ") for s in time_cols]
    size_col = "Size"
    df[time_cols] /= 1000.0
    df[size_col] = np.rint(df[size_col]).astype(int)
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    df = df.sort_values("Size")
    axf = ax.flatten()
    axindex = 0
    axf[axindex].plot(df["Load"], df[size_col], marker="o", markersize=12, linewidth=2, label="Load File")
    axf[axindex].set_ylabel("File Size (MB)")
    axindex += 1
    axf[axindex].plot(df["Index"], df[size_col], marker="+", markersize=12, linewidth=2, label="Index File")
    axf[axindex].set_ylabel("File Size (MB)")
    # df = df.sort_values('N_chan')
    # axf[2].plot(df['Create_Obsblocks'],df['N_chan'],
    #    marker='>',markersize=12,linewidth=2,label='Create Obsblocks')
    df["multi"] = df["N_chan"] * df["N_rows"] / 1e8
    df = df.sort_values("multi")
    axindex += 1
    l1 = axf[axindex].plot(
        df["Create_Obsblocks"], df["multi"], marker="^", markersize=12, linewidth=2, label="Create Obsblocks (scaled)"
    )
    axf[axindex].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    ax2 = axf[axindex].twinx()
    ax2.set_ylabel(r"$N_{rows}$")
    l2 = ax2.plot(
        df["Create_Obsblocks"],
        df["N_rows"],
        marker="o",
        markersize=12,
        linewidth=2,
        color="orange",
        label="Create Obsblocks (row)",
    )
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    axsecond = axindex
    df = df.sort_values("N_chan")
    axindex += 1
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(
        df["Create_Obsblocks"] / df["N_rows"],
        df["N_chan"],
        marker="<",
        markersize=12,
        linewidth=2,
        label="Create Obsblocks per spectrum",
    )
    axindex += 1
    axf[axindex].plot(
        df["Baseline_1"] / df["N_rows"],
        df["N_chan"],
        marker=">",
        markersize=12,
        linewidth=2,
        label="Baseline 1 per spectrum",
    )
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(
        df["Baseline_2"] / df["N_rows"],
        df["N_chan"],
        marker="+",
        markersize=12,
        linewidth=2,
        label="Baseline 2 per spectrum",
    )
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(
        df["Baseline_3"] / df["N_rows"],
        df["N_chan"],
        marker="o",
        markersize=12,
        linewidth=2,
        label="Baseline 3 per spectrum",
    )
    axf[axindex].set_ylabel(r"$N_{chans}$")
    # axf[0].set_xlabel("Time (s)")
    # axf[1].set_xlabel("Time (s)")
    # axf[2].set_xlabel("Time (s)")
    #    ax[0].set_xlabel("Time (s)")
    # axf[2].set_ylabel(r"$N_{chans}$")
    # axf[3].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    # axf[4].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    # axf[5].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    ax2.ticklabel_format(axis="y", style="sci", useMathText=True, scilimits=(0, 0))
    axindex += 1
    for j in range(axindex):
        axf[j].set_xlabel("Elapsed Time (s)")
        axf[j].ticklabel_format(axis="y", style="sci", useMathText=True, scilimits=(0, 0))
        if j == axsecond:
            axf[j].legend(lines, labels, loc="upper left")
        else:
            axf[j].legend(loc="upper left")
        # j.set_xscale('log')
        # j.set_yscale('log')
    for j in range(axindex, len(axf)):
        axf[j].axis("off")
    plt.subplots_adjust(wspace=0.35, hspace=0.25)
    fontdict = {"size": 14, "fontweight": "bold"}
    fig.suptitle(args.title, size=14, weight="bold")
    if args.outfile:
        plt.savefig(args.outfile, dpi=300)
    plt.show()


def barplots(file, x_col="N_rows"):
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.minor.size"] = 4
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.minor.size"] = 4
    plt.rcParams["font.size"] = 22
    plt.rcParams["axes.linewidth"] = 1.5
    t_python = Table.read(file[0], format="ipac")
    t_idl = Table.read(file[1], format="ipac")
    colors = [
        "#377eb8",
        "#ff7f00",
        "#4daf4a",
        "#f781bf",
        "#a65628",
        "#984ea3",
        "#999999",
        "#e41a1c",
        "#dede00",
        "#595959",
        "#5F9ED1",
        "#C85200",
        "#898989",
        "#A2C8EC",
        "#FFBC79",
        "#CFCFCF",
    ]

    size_col = "Size"
    df = dict()
    df["p"] = t_python.to_pandas().sort_values(size_col)
    df["i"] = t_idl.to_pandas().sort_values(size_col)
    time_cols = ["Load", "Index", "Create_Obsblocks", "Baseline_1", "Baseline_2", "Baseline_3"]
    tc = [s.replace("_", "\n") for s in time_cols]
    tc[2] = "Create\nSpectra"
    # df[time_cols] /= 1000.0
    for key in df:
        df[key][size_col] = np.rint(df[key][size_col]).astype(int)
        for j in time_cols:
            df[key][j] = df[key][j] / df[key]["N_rows"]

    fig, ax = plt.subplots(figsize=(12, 10))
    # print("TIME ",list(df["p"][time_cols].mean()))
    # print(df["p"])
    ind = np.arange(len(time_cols))
    width = 0.2
    ax.bar(ind, list(df["p"][time_cols].mean()), width, log=args.logy, label="dysh prototype", color=colors[6])
    ax.bar(ind + width, list(df["i"][time_cols].mean()), width, log=args.logy, label="GBTIDL", color=colors[4])
    ax.set_xticks(ind + width, labels=tc)
    ax.set_ylabel("Operation Time per spectrum (ms)", fontweight="demibold", labelpad=6)
    ax.set_xlabel("Operation", fontweight="demibold", labelpad=6)
    pure_numpy = 2.5
    pure_numpy = {0: 0.0, 1: 0, 2: 0, 3: 3.6, 4: 7.2, 5: 10.86}
    ax.bar(ind + 2 * width, list(pure_numpy.values()), width, log=args.logy, label="pure numpy", color=colors[0])
    # from peter
    pure_c = {0: 0.6, 1: 0, 2: 0, 3: 0.58, 4: 0.81, 5: 1.16}
    pure_c = {0: 0.0, 1: 0, 2: 0, 3: 0.58, 4: 0.81, 5: 1.16}
    # take 0.75 because peter did the whole spectrum
    ax.bar(
        ind + 3 * width, 0.75 * np.array(list(pure_c.values())), width, log=args.logy, label="pure C", color=colors[8]
    )

    # ax.hlines(pure_numpy,2.75,4.,lw=3,color=colors[0],label="pure numpy")
    handles, labels = ax.get_legend_handles_labels()
    # print(handles)
    # print(sorted(labels))
    # sort both labels and handles by labels
    # print(sorted(zip(labels,handles),key=lambda t: t[0].lower()))
    # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    #    ax.legend([handles[1],handles[2],handles[3],handles[0]], [labels[1],labels[2],labels[3],labels[0]])
    ax.legend()
    ax.set_title(args.title, weight="bold")
    plt.tight_layout()
    if args.show:
        plt.show()
    if args.outfile:
        plt.savefig(args.outfile, dpi=150)


# --------------------------------------------
# ax = df["p"].plot.bar(x=size_col, y=time_cols,
#        stacked=False, logy=False,
#        xlabel='File Size', ylabel = 'Time per spectrum (ms)',rot=0)#figsize=(15,12),rot=0)
# fig = ax.get_figure()
# ax2 = df["i"].plot.bar(x=size_col, y=time_cols,
#        stacked=False, logy=False, ax=ax,rot=0)#figsize=(15,12),rot=0)
#    print(ax.containers)
#    ax.bar_label(ax.containers[2],labels=df[size_col],backgroundcolor='blue',color='white',padding=-30)
#    ax.bar_label(ax.containers[2],labels=df['N_chan'],backgroundcolor='gray',color='white',padding=-50)
#    legend = ax.get_legend()
#    #legend = add_patch(legend,fc='blue',label='File Size (MB)')
#    legend = add_patch(legend,fc=['blue','gray'],
#             label=['File Size (MB)','# Channels'])
# ax.set_title("Timing for SDFITSLoad operations")

# plt.savefig("sdfitsload_timing.png",dpi=300)
#    ax = df.plot.bar(x=x_col, y=time_cols[0:2],
#            stacked=False, logy=False,
#            xlabel='Number of Rows', ylabel = 'Time (s)',rot=0)#figsize=(15,12),rot=0)
#    ax.bar_label(ax.containers[-1],labels=df[size_col],backgroundcolor='blue',color='white')
#    ax.bar_label(ax.containers[-1],labels=df['N_chan'],backgroundcolor='gray',color='white',padding=20)
#    legend = ax.get_legend()
# legend = add_patch(legend,fc='blue',label='File Size (MB)')
#    legend = add_patch(legend,fc=['blue','gray'],
#             label=['File Size (MB)','# Channels'])
#    #plt.subplot_tool(targetfig=fig)
# plt.savefig("sdfitsload_timing_short.png",dpi=300)

# ax.tick_params(axis='x',rotation=0)
# ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="revisedstructure")
    parser.add_argument("--file", "-f", action="store", help="input filename")
    parser.add_argument("--file2", "-g", action="store", help="2nd input filename [gbtidl], barplots only")
    parser.add_argument("--barplots", "-b", action="store_true", help="show barplots", default=False)
    parser.add_argument("--lineplots", "-l", action="store_true", help="show lineplots", default=False)
    parser.add_argument(
        "--title", "-t", action="store", help="Plot title", default="Timing for Load/Obsblocks/Baseline"
    )
    parser.add_argument("--outfile", "-o", action="store", help="output file", default=None)
    parser.add_argument("--logy", action="store_true", help="log y axis", default=False)
    parser.add_argument("--logx", action="store_true", help="log x axis", default=False)
    parser.add_argument("--show", action="store_true", help="show plot", default=False)
    args = parser.parse_args()
    print("ARGS ", args)
    if args.lineplots:
        lineplots(file=args.file)
    if args.barplots:
        barplots(file=[args.file, args.file2])
if False:
    for c in ax.containers:
        # Optional: if the segment is small or 0, customize the labels
        labels = [np.rint(v.get_height()).astype(int) if v.get_height() > 0 else "" for v in c]

        # remove the labels parameter if it's not needed for customized labels
        ax.bar_label(c, labels=labels, label_type="center")
