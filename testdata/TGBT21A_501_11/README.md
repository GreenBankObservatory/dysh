# TGBT21A_501_11

Test data from TGBT21A_501_11

## Files

In this directory you should find:

* NGC2782/
* NGC2782_blanks/
* testselection.fits - used to test Selection
* TGBT21A_501_11_getps_scan_152_ifnum_0_plnum_0_smthoff_15.fits
* TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits
* TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits
* TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0.fits
* TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits
* TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_0.fits 
     - used to test gettp
* TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_1.fits 
    - used to test gettp
* TGBT21A_501_11_ifnum_0_int_0-2.fits
* TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits
* TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_1.fits
* TGBT21A_501_11.raw.156.fits
* TGBT21A_501_11.raw.vegas.fits
    - used to test getps and gettp
* TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits
    -  used to test getps and gettp

## Generating the test data
To generate the calibrated total power data used for comparison:
``` IDL
filein,"TGBT21A_501_01.raw.vegas"
fileout,"aTGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
gettp,152,ifnum=0,plnum=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"aTGBT21A_501_11_scan_152_ifnum_0_plnum_0_cal_state_0.fits"
gettp,152,ifnum=0,plnum=1,cal_state=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"aTGBT21A_501_11_scan_152_ifnum_0_plnum_0_cal_state_1.fits"
print,!g.s[0].tsys,!g.s[0].exposure
gettp,152,ifnum=0,plnum=1,cal_state=1
keep
```
