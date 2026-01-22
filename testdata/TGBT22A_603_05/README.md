# TGBT22A_603_05

Argus frequency switched observations.

## Making the test data

Used GBTIDL to generate the input file.

To generate the input file:

```IDL
filein,"/home/scratch/psalas/support/UMD-py/test_data/TGBT22A_603_05.raw.vegas"
fileout,"TGBT22A_603_05.raw.vegas/TGBT22A_603_05.raw.vegas.testtrim.fits"
fdnum=2
gettp,10,fdnum=fdnum
keep
gettp,11,fdnum=fdnum
keep
gettp,12,fdnum=fdnum,sig=1
keep
gettp,12,fdnum=fdnum,sig=0
keep
```
