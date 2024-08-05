# from dysh.plot.iPythonPlotter import *
import os
import threading

from dysh.plot.canvas import SpectrumPlot
from dysh.plot.renderer import EnvironmentInfo


def plot(spectrum=None):
    env_info = EnvironmentInfo()
    env_info.info()
    renderer = env_info.get_environment()

    if renderer == "Python Script":
        print("yep, Python Script")
        my_plot = SpectrumPlot(spectrum)
        my_plot.show()
    elif renderer == "Jupyter Notebook":
        print("yep, Jupyter Notebook")
    elif renderer == "IPython Shell":
        print("yep, IPython Shell")
    else:
        print("Unknown Python configuration")


if __name__ == "__main__":
    from pathlib import Path

    from dysh.fits.gbtfitsload import GBTFITSLoad

    filename = Path(
        "/home/sandboxes/vcatlett/repos/github/vcatlett/dysh/testdata/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
    )
    sdfits = GBTFITSLoad(filename)
    psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
    ta = psscan.timeaverage(weights="tsys")

    plot(ta)
