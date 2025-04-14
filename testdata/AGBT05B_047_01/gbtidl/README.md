# GBTIDL outputs for AGBT05B_047_01

To generate the output, in `GBTIDL`:

```IDL
filein,'AGBT05B_047_01.raw.acs.fits'
getps,51
fileout,'AGBT05B_047_01.getps.acs.fits'
keep
```

For getsigref tests:
```

filein,"AGBT05B_047_01.raw.acs.fits"
getsigref,53,52,ifnum=0,plnum=0,fdnum=0,eqweight=1
fileout,"getsigref_53_52_eqweight.fits"
keep
getsigref,53,52,ifnum=0,plnum=0,fdnum=0,eqweight=0
fileout,"getsigref_53_52_tsysweight.fits"
keep
```
