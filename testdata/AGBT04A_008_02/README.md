# AGBT04A_008_02

HI-survey observations.

This data contains an SDFITS with multiple binary tables.


## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/hi-survey/data/AGBT04A_008_02.raw.acs"
fileout,"AGBT04A_008_02.raw.acs/AGBT04A_008_02.raw.acs.testrim.fits"
gettp,220,cal=1
keep
gettp,220,cal=0
keep
gettp,221,cal=1
keep
gettp,221,cal=0
keep
gettp,263,cal=1
keep
gettp,263,cal=0
keep
gettp,264,cal=1
keep
gettp,264,cal=0
keep
gettp,274,cal=1
keep
gettp,274,cal=0
keep
gettp,269,cal=1
keep
gettp,269,cal=0
keep
gettp,295
keep
gettp,296
keep
```

To produce the GBTIDL results:
```IDL
filein,"AGBT04A_008_02.raw.acs"
fileout,"AGBT04A_008_02.cal.acs.testtrim.fits"
getps,220
keep
getps,263
keep
getsigref,296,295,tsys=28
keep
getsigref,296,295
keep
```
