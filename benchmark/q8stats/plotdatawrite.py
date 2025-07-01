#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

t = Table.read("datawrite.tab", format="ascii.ecsv")
mask = ~((t["name"] == "load") | (t["name"] == "getps"))
num = []
msavg = []
msmed = []
perscan = []
for i in range(1, 14):
    m2 = t["nwrite"] == i
    num.append(i)
    msavg.append(np.mean(t[mask & m2]["time"]))
    msmed.append(np.median(t[mask & m2]["time"]))

num = np.array(num)
msavg = np.array(msavg)
msmed = np.array(msmed)
z = np.polyfit(num, msmed, deg=1)
print(f"fit={z} zeropoint={np.poly1d(z)(0)}")
perscan = (msmed - z[1]) / num
print(f"average per scan = {np.mean(perscan)}")
fig, ax = plt.subplots()
ax.plot(num, msavg, color="blue", marker="o", linestyle="solid", linewidth=2, markersize=10, label="Average")
ax.plot(num, msmed, color="green", marker="^", linestyle="solid", linewidth=2, markersize=10, label="Median")
# actually need to fit and subtract off zero point "setup" time
ax.plot(num, perscan, color="black", linestyle="dashed", linewidth=2, label="Median/scan minus fixed setup time")
ax.set_ylabel("Write Time (ms)")
ax.set_xlabel("Number of Scans in Scanblock")
ax.set_title("dysh Scanblock write performance")
ax.legend(loc="best")
rcParams = {}
rcParams["xtick.major.size"] = 7
rcParams["xtick.minor.size"] = 4
rcParams["ytick.major.size"] = 7
rcParams["ytick.minor.size"] = 4
rcParams["font.size"] = 12
rcParams["axes.linewidth"] = 1.5
plt.rcParams.update(rcParams)
plt.show()
fig.savefig("datawrite.png", dpi=300)
