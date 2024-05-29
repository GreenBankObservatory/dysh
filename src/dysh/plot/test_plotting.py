from dysh.fits.gbtfitsload import GBTFITSLoad
import astropy.units as u

filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
sdfits = GBTFITSLoad(filename)
psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
ta = psscan.timeaverage(weights='tsys')
ta.plot()