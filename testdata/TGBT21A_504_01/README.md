# TGBT21A_504_01

Test data from TGBT21A_504_01: GBT Quick Reference Guide Examples.

## Files

In this directory you should find:

* TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits   - 88 rows: 1 scan, 11 int, 2 pol, 2 sig, 2 cal
* TGBT21A_504_01.cal.vegas.fits - 2 rows:   2 pol
* TGBT21A_504_01.nofold.vegas.fits - 2 rows:  2 pol

## Generating the test data

To save space in the repo, we use only one of the scans that contains frequency switched observations (20).
Note in the raw data plnum=1 is actually the first polarization (rows 0..3) and plnum=0 the 2nd set of four
(rows 4..7).

``` bash
sdfits -scans=20 -backends=vegas TGBT21A_504_01/ScanLog.fits
```

To generate the calibrated folded data used from comparison, in both polarizations:
``` IDL
filein,"TGBT21A_504_01.raw.vegas"
getfs,20,ifnum=0,plnum=0
fileout,"TGBT21A_504_01.cal.vegas.fits"
keep
getfs,20,ifnum=0,plnum=1
keep
```
To generate the calibrated UNfolded data used from comparison, in both polarizations:
``` IDL
filein,"TGBT21A_504_01.raw.vegas"
getfs,20,ifnum=0,plnum=0,/nofold
fileout,"TGBT21A_504_01.nofold.vegas.fits"
keep
getfs,20,ifnum=0,plnum=1,/nofold
keep
```

Using intnum=0 one can also test just a single  (first in this case) integration. 
GBTIDL cannot (easily?) get a select number, it's one or all.

See also notebooks/developer/proto_getfs.ipynb   and proto_getfs.py
