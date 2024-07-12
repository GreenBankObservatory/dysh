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
To generate the calibrated total power data used for comparison. These data have only sig=T (sig_state=1):

``` IDL
filein,"TGBT21A_501_01.raw.vegas"
fileout,"TGBT21A_501_11_gettp+scan_152_ifnum_0_plnum_0.fits"
gettp,152,ifnum=0,plnum=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_0.fits"
gettp,152,ifnum=0,plnum=0,cal_state=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_1.fits"
print,!g.s[0].tsys,!g.s[0].exposure
gettp,152,ifnum=0,plnum=0,cal_state=1
keep
```

For the sig=F (sig_state=0) data, switch to project TGBT21A_504_01

``` IDL
filein,"TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
fileout,"TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_1.fits"
gettp,20,ifnum=0,plnum=1,sig_state=0,cal_state=1
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_0.fits"
gettp,20,ifnum=0,plnum=1,sig_state=0,cal_state=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
fileout,"TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_all.fits"
gettp,20,ifnum=0,plnum=1,sig_state=0
print,!g.s[0].tsys,!g.s[0].exposure
keep
```

Tsys and exposure values for sig=F, scan=20, ifnum=0, plnum=1

      24.9532       14.260588
      24.9532       13.751685
      24.9532       28.012272
