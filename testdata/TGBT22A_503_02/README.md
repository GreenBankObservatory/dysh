# TGBT22A_503_02


## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/nod-KFPA/data/TGBT22A_503_02.raw.vegas"  
fileout,"TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.testtrim.fits"
gettp,62,ifnum=0,plnum=0,fdnum=2,cal=1
keep
gettp,62,ifnum=0,plnum=0,fdnum=6,cal=1
keep
gettp,62,ifnum=0,plnum=0,fdnum=2,cal=0
keep
gettp,62,ifnum=0,plnum=0,fdnum=6,cal=0
keep
gettp,63,ifnum=0,plnum=0,fdnum=2,cal=1
keep
gettp,63,ifnum=0,plnum=0,fdnum=6,cal=1
keep
gettp,63,ifnum=0,plnum=0,fdnum=2,cal=0
keep
gettp,63,ifnum=0,plnum=0,fdnum=6,cal=0
keep
```

Output file with the calibrated data:

```IDL
filein,"TGBT22A_503_02.raw.vegas.testtrim.fits"
fileout,"TGBT22A_503_02.cal.vegas.fits"
getsigref,62,63
keep
getsigref,63,62,fdnum=6
keep
```
