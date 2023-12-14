# GBTIDL spectra for testing

All of the raw data used to generate these examples should be available in the [dysh data server](https://www.gb.nrao.edu/dysh/example_data/).

Here we provide the `GBTIDL` commands used to generate individual files.


## OnOff L <a name="onoff"></a>
``` IDL
filein,"onoff-L/data/TGBT21A_501_11.raw.vegas/TGBT21A_501_11.raw.vegas.A.fits"
gettp,156,intnum=0
setframe,"TOPO"
write_ascii,"onoff-L_gettp_156_intnum_0_TOPO.ascii"
setframe,"GEO"
write_ascii,"onoff-L_gettp_156_intnum_0_GEO.ascii"
setframe,"HEL"
write_ascii,"onoff-L_gettp_156_intnum_0_HEL.ascii"
setframe,"BAR"
write_ascii,"onoff-L_gettp_156_intnum_0_BAR.ascii"
setframe,"LSR"
write_ascii,"onoff-L_gettp_156_intnum_0_LSR.ascii"
setframe,"LSD"
write_ascii,"onoff-L_gettp_156_intnum_0_LSD.ascii"
setframe,"GAL"
write_ascii,"onoff-L_gettp_156_intnum_0_GAL.ascii"
```
