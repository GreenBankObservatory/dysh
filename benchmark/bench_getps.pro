;    ln -s $DATA_DYSH/example_data/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits
filein,'AGBT05B_047_01.raw.acs.fits'
getps,999
freeze
TIC
getps,51
getps,53
getps,55
getps,57
aver
dt = TOC()
PRINT, 'Elapsed time in sec: ',dt
EXIT
