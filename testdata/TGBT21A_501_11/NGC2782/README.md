Test data including two OnOff scans.

To generate the SDFITS files, in `Python`

```
import numpy as np

from astropy.io import fits

from dysh import util
from dysh.fits.gbtfitsload import GBTFITSLoad


# Select the VEGAS banks to use,
# and the scan numbers you are interested in.
vbanks = ["A"]
scans  = [156, 157, 158, 159]

# Where do you find the original data?
example_data = "/home/dysh/example_data/"     # @GBO
data_path    = f"{example_data}onoff-L/data/TGBT21A_501_11.raw.vegas/"
data_file    = lambda vbank : f"{data_path}/TGBT21A_501_11.raw.vegas.{vbank}.fits"
# Where will you put the smaller data?
# Use meaningful names.
output_path  = f"{util.get_project_testdata()}/TGBT21A_501_11/NGC2782/"
output       = lambda vbank : f"{output_path}/TGBT21A_501_11_NGC2782.raw.vegas.{vbank}.fits"


orows  = []

# Treat each VEGAS bank independently.
for vbank in vbanks:

    # Define the input and output.
    df = data_file(vbank)
    of = output(vbank)

    # Load the input.
    sdf = GBTFITSLoad(df)

    # We want to keep at least one noise diode cycle
    # for each spectral window and polarization.
    nrows  = len(sdf.udata("PLNUM")) * len(sdf.udata("IFNUM")) * len(sdf.udata("CAL"))

    # Loop over scans getting the row number
    # of the first element.
    # @TODO: check if plnum=1 is always first.
    for scan in scans:
        rows = sdf.scan_rows([scan], plnum=1)
        print(f"First row of scan {scan} is {rows[0]}")
        orows.append(np.arange(rows[0], rows[0]+nrows))

    # Save the rows for one VEGAS bank.
    hdu0  = sdf._sdf[0]._hdu[0].copy()
    table = sdf._sdf[0]._hdu[1].data[np.ravel(orows)]
    head  = sdf._sdf[0]._hdu[1].header
    thdu  = fits.BinTableHDU(table, header=head)
    outhdu = fits.HDUList([hdu0, thdu])
    outhdu.writeto(of, overwrite=True)
    print(f"Saved rows to {of}")
```

To generate the `GBTIDL` output, working from
testdata/TGBT21A_501_11

```
filein,"NGC2782"
getps,156,ifnum=0,plnum=0
accum
getps,158,ifnum=0,plnum=0
accum
ave
fileout,"NGC2782/TGBT21A_501_11_getps_scans_156-158_ifnum_0_plnum_0_timeaverage.fits"
keep
```
