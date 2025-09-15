# TRFI_090125_S1

RFI scan observations.

This data contains an SDFITS with multiple binary tables and repeated scan numbers.

## Making the test data

Used GBTIDL to generate the input file.

To generate the input file:

```IDL
offline,'TRFI_090125_S1'
fileout,'TRFI_090125_S1.raw.vegas/TRFI_090125_S1.raw.vegas.testtrim.fits'
a=getchunk(scan=2, ifnum=0, plnum=0, int=[1,2,3,4])
putchunk,a
```
