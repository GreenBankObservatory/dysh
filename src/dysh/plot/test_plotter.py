from dysh.fits.gbtfitsload import GBTFITSLoad

# from dysh.plot.specplot import SpectrumPlot

filename = "E:/Code/GitHub/Work/forks/dysh/test_file.fits"
sdfits = GBTFITSLoad(filename)
print(sdfits.summary(show_index=True))

psscan = sdfits.getps(152, ifnum=0, plnum=0)
ta = psscan.timeaverage(weights="tsys")
# breakpoint()
ta.plot(title="Sample Plot", xaxis_unit="km/s", yaxis_unit="mK", ymin=-100, ymax=500, xmin=3000, xmax=4500)

# [TODO]
# See this for ideas: https://pyspeckit.readthedocs.io/en/latest/interactive.html
# Edit selection object with clicks
# Have file info listed at the top
# Make sure it works in scripts, iPython, and Jupyter
