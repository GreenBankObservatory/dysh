import socket

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.plot.dev.iPlotter import SpectrumPlot as sp

# Get path to data
hostname = socket.gethostname()
if hostname == "VCatlett-Desktop":
    # (Cat) For testing at home. Will remove when prod ready
    filename = "E:/Data/GBO/dysh/test_file.fits"
else:
    filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"

# Load data
sdfits = GBTFITSLoad(filename)
psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
ta = psscan.timeaverage(weights="tsys")

# Make sure .plot() works
myplot = ta.plot()

# Make sure selection works
my_selection = sp(ta)

print(my_selection.get_selection())
