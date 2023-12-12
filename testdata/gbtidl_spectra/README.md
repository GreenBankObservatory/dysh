# GBTIDL spectra for testing

All of the raw data used to generate these examples should be available in the [dysh data server](https://www.gb.nrao.edu/dysh/example_data/).

Here we provide the commands used to generate individual files.

## Table of contents
1. [spider-C](#spider)
2. [Some paragraph](#paragraph1)
    1. [Sub paragraph](#subparagraph1)
    3. [Another paragraph](#paragraph2)


## spider-C <a name="spider"></a>
``` IDL
filein,"spider-C/data/AGBT20B_424_02.raw.vegas/AGBT20B_424_02.raw.vegas.A.fits"
gettp,136,intnum=0
write_ascii,"spider-C_gettp_136_intnum_0.ascii"
```
