# AGBT21B_024_01

## Making the test data

Used GBTIDL to generate the input and output files.

To generate the input file:

```IDL
filein,"/home/dysh/acceptance_testing/data/AGBT21B_024_01/AGBT21B_024_01.raw.vegas"
fileout,"AGBT21B_024_01.raw.vegas/AGBT21B_024_01.raw.vegas.testtrim.fits"
gettp,19,fdnum=10
keep
gettp,20,fdnum=10
keep
gettp,104,fdnum=10
keep
gettp,105,fdnum=10
keep
```
