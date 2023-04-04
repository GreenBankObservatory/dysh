#!/usr/bin/env python
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import argparse

def add_patch(legend,fc,label):
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
    t = Table.read(file,format='ipac')
    colors = ['red', 'tan', 'lime']
    df = t.to_pandas()#.sort_values('N_rows')
    time_cols = ['Load','Create_Obsblocks' ,'Baseline_1','Baseline_2','Baseline_3']
    size_col = 'Size'
    df[time_cols] /= 1000.0
    df[size_col] = np.rint(df[size_col]).astype(int)
    fig,ax = plt.subplots(6,1)
    df = df.sort_values('Load')
    axf = ax.flatten()
    axf[0].plot(df['Load'],df[size_col],
        marker='o',markersize=12,linewidth=2,label='Load File')
    df = df.sort_values('N_rows')
    #axf[1].plot(df['Create_Obsblocks'],df['N_rows'],
    #    marker='<',markersize=12,linewidth=2,label='Create Obsblocks')
    #df = df.sort_values('N_chan')
    #axf[2].plot(df['Create_Obsblocks'],df['N_chan'],
    #    marker='>',markersize=12,linewidth=2,label='Create Obsblocks')
    df['multi'] = df['N_chan']*df['N_rows']/1E8
    df = df.sort_values('multi')
    axf[1].plot(df['Create_Obsblocks'],df['multi'],
        marker='^',markersize=12,linewidth=2,label='Create Obsblocks')
    axf[2].plot(df['Baseline_1'],df['multi'],
        marker='>',markersize=12,linewidth=2,label='Baseline 1')
    axf[3].plot(df['Baseline_2'],df['multi'],
        marker='>',markersize=12,linewidth=2,label='Baseline 2')
    axf[4].plot(df['Baseline_3'],df['multi'],
        marker='>',markersize=12,linewidth=2,label='Baseline 3')
    axf[0].set_xlabel("Time (s)")
    axf[1].set_xlabel("Time (s)")
    axf[2].set_xlabel("Time (s)")
#    ax[0].set_xlabel("Time (s)")
    axf[0].set_ylabel("File Size (MB)")
    #ax[1].set_ylabel(r"$N_{rows}$")
    #ax[1][0].set_ylabel(r"$N_{chans}$")
    axf[1].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    #ax[2][0].set_ylabel(r"$N_{chans}$")
    axf[2].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    #ax[2][1].set_ylabel(r"$N_{chans}$")
    axf[3].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    #ax[2][2].set_ylabel(r"$N_{chans}$")
    axf[2].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    for j in axf:
        j.legend()
        #j.set_xscale('log')
        #j.set_yscale('log')
    plt.show()
    

def barplots(file,x_col='N_rows'):
    t = Table.read(file,format='ipac')
    colors = ['red', 'tan', 'lime']
    df = t.to_pandas().sort_values('N_rows')
    time_cols = ['Load','Create_Obsblocks' ,'Baseline_1','Baseline_2','Baseline_3']
    size_col = 'Size'
    df[time_cols] /= 1000.0
    df[size_col] = np.rint(df[size_col]).astype(int)

    ax = df.plot.bar(x=x_col, y=time_cols,
            stacked=False, logy=False, 
            xlabel='Number of Rows', ylabel = 'Time (s)',rot=0)#figsize=(15,12),rot=0)
    print(ax.containers)
    ax.bar_label(ax.containers[2],labels=df[size_col],backgroundcolor='blue',color='white',padding=-30)
    ax.bar_label(ax.containers[2],labels=df['N_chan'],backgroundcolor='gray',color='white',padding=-50)
    legend = ax.get_legend()
    #legend = add_patch(legend,fc='blue',label='File Size (MB)')
    legend = add_patch(legend,fc=['blue','gray'],
             label=['File Size (MB)','# Channels'])
    ax.set_title("Timing for SDFITSLoad operations")

    #plt.savefig("sdfitsload_timing.png",dpi=300)
    plt.show()
    ax = df.plot.bar(x=x_col, y=time_cols[0:2],
            stacked=False, logy=False, 
            xlabel='Number of Rows', ylabel = 'Time (s)',rot=0)#figsize=(15,12),rot=0)
    ax.bar_label(ax.containers[-1],labels=df[size_col],backgroundcolor='blue',color='white')
    ax.bar_label(ax.containers[-1],labels=df['N_chan'],backgroundcolor='gray',color='white',padding=20)
    legend = ax.get_legend()
    #legend = add_patch(legend,fc='blue',label='File Size (MB)')
    legend = add_patch(legend,fc=['blue','gray'],
             label=['File Size (MB)','# Channels'])
    ax.set_title("Timing for SDFITSLoad operations")
    plt.show()
    #plt.savefig("sdfitsload_timing_short.png",dpi=300)

#ax.tick_params(axis='x',rotation=0)
#ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="revisedstructure")
    parser.add_argument('--file','-f', action='store', help='input filename')
    parser.add_argument('--barplots','-b', action='store_true', help='show barplots',default=False)
    parser.add_argument('--lineplots','-l', action='store_true', help='show lineplots',default=True)
    args = parser.parse_args()
    print("ARGS ",args)
    if args.lineplots:
        lineplots(file=args.file)
    if args.barplots:
        barplots(file=args.file)
if False:
    for c in ax.containers:

        # Optional: if the segment is small or 0, customize the labels
        labels = [np.rint(v.get_height()).astype(int) if v.get_height() > 0 else '' for v in c]

        # remove the labels parameter if it's not needed for customized labels
        ax.bar_label(c, labels=labels, label_type='center')
