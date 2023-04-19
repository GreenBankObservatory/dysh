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
    time_cols = ['Load','Index', 'Create_Obsblocks' ,'Baseline_1','Baseline_2','Baseline_3']
    size_col = 'Size'
    df[time_cols] /= 1000.0
    df[size_col] = np.rint(df[size_col]).astype(int)
    fig,ax = plt.subplots(3,2,figsize=(10,10))
    df = df.sort_values('Size')
    axf = ax.flatten()
    axindex = 0
    axf[axindex].plot(df['Load'],df[size_col],
        marker='o',markersize=12,linewidth=2,label='Load File')
    axf[axindex].set_ylabel("File Size (MB)")
    axindex+=1
    axf[axindex].plot(df['Index'],df[size_col],
        marker='+',markersize=12,linewidth=2,label='Index File')
    axf[axindex].set_ylabel("File Size (MB)")
    #df = df.sort_values('N_chan')
    #axf[2].plot(df['Create_Obsblocks'],df['N_chan'],
    #    marker='>',markersize=12,linewidth=2,label='Create Obsblocks')
    df['multi'] = df['N_chan']*df['N_rows']/1E8
    df = df.sort_values('multi')
    axindex+=1
    l1 = axf[axindex].plot(df['Create_Obsblocks'],df['multi'],
        marker='^',markersize=12,linewidth=2,
        label='Create Obsblocks (scaled)')
    axf[axindex].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    ax2 = axf[axindex].twinx()
    ax2.set_ylabel(r"$N_{rows}$")
    l2 = ax2.plot(df['Create_Obsblocks'],df['N_rows'],
        marker='o',markersize=12,linewidth=2,color='orange',
        label='Create Obsblocks (row)')
    lines = l1+l2
    labels = [l.get_label() for l in lines]
    axsecond = axindex
    df = df.sort_values('N_chan')
    axindex+=1
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(df['Create_Obsblocks']/df['N_rows'],df['N_chan'],
        marker='<',markersize=12,linewidth=2,label='Create Obsblocks per spectrum')
    axindex+=1
    axf[axindex].plot(df['Baseline_1']/df['N_rows'],df['N_chan'],
        marker='>',markersize=12,linewidth=2,label='Baseline 1 per spectrum')
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(df['Baseline_2']/df['N_rows'],df['N_chan'],
        marker='+',markersize=12,linewidth=2,label='Baseline 2 per spectrum')
    axf[axindex].set_ylabel(r"$N_{chans}$")
    axf[axindex].plot(df['Baseline_3']/df['N_rows'],df['N_chan'],
        marker='o',markersize=12,linewidth=2,label='Baseline 3 per spectrum')
    axf[axindex].set_ylabel(r"$N_{chans}$")
    #axf[0].set_xlabel("Time (s)")
    #axf[1].set_xlabel("Time (s)")
    #axf[2].set_xlabel("Time (s)")
#    ax[0].set_xlabel("Time (s)")
    #axf[2].set_ylabel(r"$N_{chans}$")
    #axf[3].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    #axf[4].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    #axf[5].set_ylabel(r"$N_{chans} \times * N_{rows}$ (scaled)")
    ax2.ticklabel_format(axis='y',style='sci',useMathText=True,scilimits=(0,0))
    axindex+=1
    for j in range(axindex):
        axf[j].set_xlabel("Elapsed Time (s)")
        axf[j].ticklabel_format(axis='y',style='sci',useMathText=True,scilimits=(0,0))
        if j == axsecond:
            axf[j].legend(lines,labels,loc="upper left")
        else:
            axf[j].legend(loc="upper left")
        #j.set_xscale('log')
        #j.set_yscale('log')
    for j in range(axindex,len(axf)):
        axf[j].axis('off')
    plt.subplots_adjust(wspace=0.35,hspace=0.25)
    fontdict = {'size':14,'fontweight':'bold'}
    fig.suptitle(args.title,size=14,weight='bold')
    if args.outfile:
        plt.savefig(args.outfile,dpi=300)
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
    plt.subplot_tool(targetfig=fig)
    #plt.show()
    #plt.savefig("sdfitsload_timing_short.png",dpi=300)

#ax.tick_params(axis='x',rotation=0)
#ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="revisedstructure")
    parser.add_argument('--file','-f', action='store', help='input filename')
    parser.add_argument('--barplots','-b', action='store_true', help='show barplots',default=False)
    parser.add_argument('--lineplots','-l', action='store_true', help='show lineplots',default=True)
    parser.add_argument('--title','-t', action='store', help='Plot title',default="Timing for Load/Obsblocks/Baseline")
    parser.add_argument('--outfile','-o', action='store', help='output file',default=None)
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
