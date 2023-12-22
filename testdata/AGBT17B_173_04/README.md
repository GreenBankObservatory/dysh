# AGBT17B_173_04

Test data from AGBT17B_173_04: The GBT Diffuse Ionized Gas Survey (GDIGS)

## Generating the test data

To save space in the repo, we use only the scans that contain poistion switched observations (6 & 7), select only three "good" spectral windows and average together the integrations.

Create an SDFITS file with only scans 6 and 7
``` bash
sdfits -scans=6,7 AGBT17B_173_04/ScanLog.fits
```

Average together all the integrations, in `GBTIDL`
``` IDL
filein,"AGBT17B_173_04.raw.vegas"
ifnums=[0,19,42]
fileout,"gdigs-testdata.fits"
for scan=6,7,1 do begin & for plnum=0,1,1 do begin & for i=0,2,1 do begin & for cal=0,1,1 do begin & gettp,scan,ifnum=ifnums[i],plnum=plnum,cal_state=cal & keep & endfor & endfor & endfor & endfor
```

Then, we can generate the output for the tests
``` IDL
fileout,"gdigs-testdata-getps-outputs.fits"
filein,"gdigs-testdata.fits"
for plnum=0,1,1 do begin & for i=0,2,1 do begin & getps,6,ifnum=ifnums[i],plnum=plnum & keep & velo & write_ascii,"ifnum_"+strcompress(string(ifnums[i]), /remove_all)+"_plnum_"+strcompress(string(plnum),/remove_all)+".ascii" & endfor & endfor
```
