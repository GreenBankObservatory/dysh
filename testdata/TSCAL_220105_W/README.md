# TSCAL_220105_W

## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/rxco-W/data/TSCAL_220105_W.raw.vegas"
fileout,"TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.fits"
gettp,24,plnum=0,fdnum=0
keep
gettp,24,plnum=0,fdnum=1
keep
gettp,25,plnum=0,fdnum=0
keep
gettp,25,plnum=0,fdnum=1
```

Output file with the calibrated data:

```IDL
filein,"TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.fits"
fileout,"TSCAL_220105_W.cal.vegas.fits"
getsigref,24,25,fdnum=0
keep
getsigref,25,24,fdnum=0
keep
```
