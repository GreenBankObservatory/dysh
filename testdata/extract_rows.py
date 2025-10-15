#! /usr/bin/env python3
#
# A script to extract specific rows from an SDFITS file and write them
# out to a new SDFITS file.

import sys

import numpy as np
from astropy.io import fits

from dysh.fits.gbtfitsload import GBTFITSLoad

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} row_min row_max | row1 row2 row3 .... rowN")
    print("Select either a row range (min,max), or selected rows.")
    print("Row numbers are 0 based")
    print("Example of selecting many:  $(seq  0 15)")
    print("Output filename is always junk.fits, so please rename")
    sys.exit(0)

iname = sys.argv[1]
oname = "junk.fits"

if len(sys.argv[2:]) == 2:
    rows = list(range(int(sys.argv[2]), int(sys.argv[3]) + 1))
else:
    rows = []
    for r in sys.argv[2:]:
        rows.append(int(r))
print("Selecting rows", rows)

s2 = GBTFITSLoad(iname)
hdu0 = s2._sdf[0]._hdu[0].copy()
nrows = len(s2._sdf[0]._hdu[1].data)
table = s2._sdf[0]._hdu[1].data[np.ravel(rows)]
head = s2._sdf[0]._hdu[1].header
hdu1 = fits.BinTableHDU(table, header=head)
outhdu = fits.HDUList([hdu0, hdu1])
outhdu.writeto(oname, overwrite=True)
print("Wrote", oname)
