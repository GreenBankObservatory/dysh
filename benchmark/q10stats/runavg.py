#!/usr/bin/env python
from astropy.table import Table
import numpy as np
#q8stats = "./redoq8/loadtimes_before"
q8stats = "../q8stats/loadtimes_before"
q10stats = "loadtimes_after"

t = {}
for k in [q8stats, q10stats]:
    print(k)
    p = Table.read(k,format='ascii.basic')
    #p.pprint_all()
    t[k] = Table(names=['data', 'mean_load', 'median_load', 'mean_VMRSS', '#files',  '#flags', 'skipflags'], units=["", "ms", "ms", "MByte", "", "", ""], dtype=[str,float,float,float,int,int,str])
    #print('data mean(ms) median(ms) mean(VMRSS) #files  #flags skipflags')
    for i in np.arange(len(p),step=4):
    #    print(f"{p[i]['data']} {np.mean(p[i:i+4]['time']):.1f} {np.median(p[i:i+4]['time']):.1f} {np.mean(p[i:i+4]['VmRSS']):.1f} {p[i]['#files']} {p[i]['#flags']} {p[i]['skipflags']}")

        vals=[p[i]['data'], np.mean(p[i:i+4]['time']), np.median(p[i:i+4]['time']), np.mean(p[i:i+4]['VmRSS']), p[i]['#files'], p[i]['#flags'], p[i]['skipflags']]
        #print(vals)
        t[k].add_row(vals)

    t[k].pprint_all()

tabdiff = Table(data=[
    t[q8stats]['data'], 
    t[q8stats]['mean_load'], 
    t[q10stats]['mean_load'], 
    (t[q8stats]['mean_load']-t[q10stats]['mean_load']),
    (t[q8stats]['median_load']-t[q10stats]['median_load']),
    t[q10stats]['skipflags']],
    names=["File", "Load Time Before", "Load Time After", "Mean Diff", "Median Diff", "skipflags"],
    units=["",        "ms",               "ms",               "ms",         "ms",             ""])
#print(f"DIFF : {(t[q8stats]['mean_load']-t[q10stats]['mean_load'])} {t[q8stats]['skipflags']} {t[q10stats]['skipflags']}")
#print(f"DIFF : {(t[q8stats]['median_load']-t[q10stats]['median_load'])} {t[q8stats]['skipflags']} {t[q10stats]['skipflags']}")
for c in tabdiff.columns:
    if tabdiff[c].unit == "ms":
        tabdiff[c].format='{:.1f}'
tabdiff.pprint_all()
tabdiff.write('finaldiff.csv',format='ascii.csv',overwrite=True)
