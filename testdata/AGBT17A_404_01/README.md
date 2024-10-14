# GBT17A-404

Position-switched HI scans of a nearby galaxy from project [AGBT17A-404](https://dss.gb.nrao.edu/project/GBT17A-404) (PI D.J. Pisano). This data is used for the baseline modelling tests and features a 3rd order baseline with 2 regions of interest (3 regions excised from the baseline modelling)

Some of the integrations have a significant GPS-L3 RFI spike, requiring the excision of integrations 43-51.


To generate the GBTIDL results:
```
filein,'AGBT17A_404_01.raw.vegas.A.fits'
for i=0,42 do begin & getps,19,plnum=0,intnum=i & accum & end
for i=0,42 do begin & getps,19,plnum=1,intnum=i & accum & end
for i=52,60 do begin & getps,19,plnum=0,intnum=i & accum & end
for i=52,60 do begin & getps,19,plnum=1,intnum=i & accum & end
ave
boxcar,5,/decimate
copy,0,1
fileout,'AGBT17A_404_01_scan_19_prebaseline.fits'
keep
;precisely place fitting regions for baseline
!g.nregion=2
!g.regions[0,0] = 100
!g.regions[1,0] = 380
!g.regions[0,1] = 450
!g.regions[1,1] = 720
bshape,nfit=3,modelbuffer=2
bsubtract
fileout,'AGBT17A_404_01_scan_19_postbaseline.fits'
keep
copy,2,0
fileout,'AGBT17A_404_01_scan_19_bmodel.fits'
keep
```

To generate the dysh results:
```
import numpy as np
from dysh.fits.gbtfitsload import GBTFITSLoad
import dysh

#where to find the data?
fname_base = "[/path/to/dysh]/dysh/testdata/AGBT17A_404_01/"
pre_fname = fname_base + 'AGBT17A_404_01_scan_19_prebaseline.fits'
post_fname = fname_base + 'AGBT17A_404_01_scan_19_postbaseline.fits'
bline_fname = fname_base+ 'AGBT17A_404_01_scan_19_bmodel.fits'

#load in with dysh, get data
sdf = dysh.fits.gbtfitsload.GBTFITSLoad(pre_fname)
dysh_spectrum = sdf.getspec(0)

#set baseline parameters and fit
order = 3
excised_regions = [(0,100),(381,450),(751,820)]




```

