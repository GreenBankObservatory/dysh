# TGBT24B_615_01

Argus OnOff observations.

## Making the test data

Using ``dysh``

```Python
from dysh.fits import GBTFITSLoad
fnm = "/home/scratch/psalas/support/UMD-py/test_data/TGBT24B_615_01.raw.vegas"
sdf = GBTFITSLoad(fnm, skipflags=True, flag_vegas=True)
sdf.write("TGBT24B_615_01.raw.vegas/TGBT24B_615_01.raw.vegas.testtrim.fits", scan=[84,85,86,87], ifnum=0, plnum=0, fdnum=10)
```
