# TGBT21A_504_01

Test data from TGBT21A_504_01: GBT Quick Reference Guide Examples.

## Generating the test data

To save space in the repo, we use only one of the scans that contains frequency switched observations (20).

``` bash
sdfits -scans=20 -backends=vegas TGBT21A_504_01/ScanLog.fits
```

To generate the calibrated data used from comparison:
``` IDL
filein,"TGBT21A_504_01.raw.vegas"
getfs,20,ifnum=0,plnum=0
fileout,"TGBT21A_504_01.cal.vegas.fits"
keep
getfs,20,ifnum=0,plnum=1
keep
```
