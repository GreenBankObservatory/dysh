#!/bin/csh -f
#set path = ( . $path )
#
# BE sure to echo 1 > /proc/sys/vm/drop_caches before running this.
unsetenv OMP_NUM_THREADS 
set root = "/data/gbt/examples/"
#set root = "/lma1/mpound/GBT/examples/"
#set root = "/home/scratch/mpound/examples/"
set files = ( \
"rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.E.fits"  \
#"misc/IC1481.fits"  \
#"misc/W3OH.fits"  \
#"misc/ngc5291.fits"  \
#"onoff-L/data/TGBT21A_501_11.raw.vegas.fits" \
#"mapping-L/data/TGBT17A_506_11.raw.vegas/TGBT17A_506_11.raw.vegas.A.fits"  \
#"nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.F.fits"  \
#"mixed-fs-ps/data/AGBT16B_225_05/AGBT16B_225_05.raw.vegas/AGBT16B_225_05.raw.vegas.B.fits" \
)

set date=`date +'%Y-%m-%d'`
set outfile = "stats-wcs-meta.$date-OMPunset"
rm -rf tab.out $outfile
foreach f ($files)
    set infile=$root$f
    set size=`du -s -BM $infile`
    echo $size
    python revisedstructure.py -b -i -o -f $infile >> $outfile
end
