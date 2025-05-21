#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

t = Table.read("datawrite.tab",format='ascii.ecsv')
mask = ~((t['name'] == 'load') | (t['name'] == 'getps'))
num = []
msavg = []
msmed = []
for i in range(1,14):
    m2 = t['nwrite'] == i
    num.append(i)
    msavg.append(np.mean(t[mask & m2]['time']))
    msmed.append(np.median(t[mask & m2]['time']))
    
fig,ax = plt.subplots()
ax.plot(num,msavg,color='blue',marker='o',linestyle='solid',linewidth=2,markersize=10,label='Average')
ax.plot(num,msmed,color='green',marker='^',linestyle='solid',linewidth=2,markersize=10,label='Median')
ax.set_ylabel("Write Time (ms)")
ax.set_xlabel("Number of Scans in Scanblock")
ax.set_title("dysh Scanblock write performance")
ax.legend(loc='best')
rcParams = {}
rcParams["xtick.major.size"] = 7
rcParams["xtick.minor.size"] = 4
rcParams["ytick.major.size"] = 7
rcParams["ytick.minor.size"] = 4
rcParams['font.size'] = 12
rcParams['axes.linewidth'] =1.5
plt.rcParams.update(rcParams)
plt.show()
fig.savefig("datawrite.png",dpi=300)
