;    ln -s $DATA_DYSH/example_data/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits
filein,'AGBT05B_047_01.raw.acs.fits'
freeze
TIC
; without timeaver
; for i=51,57,2 do begin getps, i & end
; with timeaver
for i=51,57,2 do begin getps, i & accum & end
ave
dt = TOC()
PRINT, 'Elapsed time in sec: ',dt
stats,6000,12000
EXIT
