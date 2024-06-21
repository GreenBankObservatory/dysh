<<<<<<< HEAD:src/dysh/plot/dev_tests/test_plotting.py
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.plot.dev.ex_plotly import SpecPlot

=======
>>>>>>> f6d3f8953d59c3dbebf6d1ca22b994c2440478a1:src/dysh/plot/test_plotting.py
import astropy.units as u

from dysh.fits.gbtfitsload import GBTFITSLoad

filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
sdfits = GBTFITSLoad(filename)
psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
<<<<<<< HEAD:src/dysh/plot/dev_tests/test_plotting.py
ta = psscan.timeaverage(weights='tsys')
ta.mask[0:300] = True
test_plot = SpecPlot(ta)
test_plot.show()
=======
ta = psscan.timeaverage(weights="tsys")
ta.plot()
>>>>>>> f6d3f8953d59c3dbebf6d1ca22b994c2440478a1:src/dysh/plot/test_plotting.py
