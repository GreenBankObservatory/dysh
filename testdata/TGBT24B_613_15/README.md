# TGBT24B_613_15

## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/onoff-Argus/data/TGBT24B_613_15.raw.vegas"
fileout,"TGBT24B_613_15.raw.vegas/TGBT24B_613_15.raw.vegas.testtrim.fits"
gettp,23,fdnum=10
keep
gettp,24,fdnum=10
keep
```

Output file with the calibrated data:

```IDL
filein,"TGBT24B_613_15.raw.vegas/TGBT24B_613_15.raw.vegas.testtrim.fits"
fileout,"TGBT24B_613_15.cal.vegas.fits"
getsigref,23,24,fdnum=10
keep
```
