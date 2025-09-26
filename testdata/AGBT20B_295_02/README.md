# AGBT20B_295_02

Argus frequency switched observations.

## Making the test data

Used GBTIDL to generate the input file.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/fs-Argus/data/AGBT20B_295_02.raw.vegas/AGBT20B_295_02.raw.vegas.A.fits"
fileout,"AGBT20B_295_02.raw.vegas.testtrim.fits"

a=getchunk(scan=10,fdnum=10)
putchunk,a
a=getchunk(scan=11,fdnum=10)
putchunk,a
a=getchunk(scan=12,fdnum=10)
putchunk,a[0:11]
```
