Test data including one OnOff scan with blanked integrations.

To generate the SDFITS files, in `Python`

```
import numpy as np
import pandas as pd

from astropy.io import fits

from dysh import util
from dysh.fits.gbtfitsload import GBTFITSLoad

# Select the VEGAS banks to use,
# and the scan numbers you are interested in.
vbanks = ["A"]
scans  = [156, 157]
nints  = 2

# Where do you find the original data?
example_data = "/home/dysh/example_data/"     # @GBO
data_path    = f"{example_data}onoff-L/data/TGBT21A_501_11.raw.vegas/"
data_file    = lambda vbank : f"{data_path}/TGBT21A_501_11.raw.vegas.{vbank}.fits"
# Where will you put the smaller data?
# Use meaningful names.
output_path  = f"{util.get_project_testdata()}/TGBT21A_501_11/NGC2782_blanks/"
output       = lambda vbank : f"{output_path}/NGC2782.raw.vegas.{vbank}.fits"

for vbank in vbanks:

    file_in  = data_file(vbank)
    file_out = output(vbank)

    # Open.
    hdu = fits.open(file_in)
    table = hdu[1].data

    # Bookkeeping.
    nif = len(list(set(table["IFNUM"])))
    npl = len(list(set(table["PLNUM"])))
    nnd = len(list(set(table["CAL"])))
    nsw = len(list(set(table["SIG"])))

    orows = []
    nrows = []

    for i,scan in enumerate(scans):

        # Bookkeeping.
        nrows.append(nif * npl * nnd * nsw)

        # Row selection.
        mask = np.all(np.isnan(table["DATA"]), axis=1) & (table["SCAN"] == scan)
        rows = np.where(mask == True)[0]
        for j in range(nints):
            skip = j*nrows[i]
            rowi = (rows[0] - 1) - skip
            rowf = rows[-1] - skip
            orows += list(np.arange(rowi, rowf+1))

    # Check results.
    # The resulting data frame should have entries for all
    # noise diode states in all polarizations, spectral windows
    # and switching states.
    table_ = np.lib.recfunctions.drop_fields(hdu[1].data, "DATA")
    df = pd.DataFrame(table_)
    for i,nrow in enumerate(nrows*nints):
        irow = i*nrow
        frow = (i+1)*nrow
        print(df[["SCAN", "IFNUM", "PLNUM", "CAL", "SIG", "EXPOSURE"]].iloc[np.sort(np.ravel(orows))[irow:frow]])
        print("\n")
    # If the results are correct
    # save the rows for one VEGAS bank.
    hdu0  = hdu[0].copy()
    table = hdu[1].data[np.sort(np.ravel(orows))]
    head  = hdu[1].header
    thdu  = fits.BinTableHDU(table, header=head)
    outhdu = fits.HDUList([hdu0, thdu])
    outhdu.writeto(file_out, overwrite=True)
    print(f"Saved rows to {file_out}")
```
