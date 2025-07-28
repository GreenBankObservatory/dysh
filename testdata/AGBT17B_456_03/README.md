# AGBT17B_456_03

Argus sub beam nod observations.

## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/example_data/Argus-sbn/data/AGBT17B_456_03.raw.vegas/AGBT17B_456_03.raw.vegas.A.fits"
fileout,"AGBT17B_456_03.raw.vegas.testtrim.fits"

a=getchunk(scan=20,fdnum=10)
putchunk,a[0:60]
```
