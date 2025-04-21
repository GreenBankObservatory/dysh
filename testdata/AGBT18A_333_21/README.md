# AGBT18A_333_21


## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/acceptance_testing/data/AGBT18A_333_21/AGBT18A_333_21.raw.vegas"
fileout,"AGBT18A_333_21.raw.vegas/AGBT18A_333_21.raw.vegas.testrim.fits"
gettp,11,ifnum=0,plnum=0,fdnum=1,cal=1                                         
keep
gettp,11,ifnum=0,plnum=1,fdnum=0,cal=1
keep
gettp,11,ifnum=0,plnum=0,fdnum=1,cal=0
keep
gettp,11,ifnum=0,plnum=1,fdnum=0,cal=0
keep
gettp,12,ifnum=0,plnum=0,fdnum=1,cal=1                                  
keep
gettp,12,ifnum=0,plnum=1,fdnum=0,cal=1 
keep
gettp,12,ifnum=0,plnum=0,fdnum=1,cal=0 
keep
gettp,12,ifnum=0,plnum=1,fdnum=0,cal=0 
keep
gettp,197,ifnum=0,plnum=0,fdnum=0,cal=1
keep
gettp,197,ifnum=0,plnum=0,fdnum=0,cal=0
keep
gettp,198,ifnum=0,plnum=0,fdnum=0,cal=1
keep
gettp,198,ifnum=0,plnum=0,fdnum=0,cal=0
keep
```
