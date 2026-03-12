#!/bin/bash
#
#   takes ~10 mins on lma,   ~7 mins on peter's laptop
#
#   Don't forget:     sync;sync;sync; echo 1 | sudo tee /proc/sys/vm/drop_caches

# Script to run various benchmarks
export OMP_NUM_THREADS=1

#   meant to run on "lma", you need to set DYSH_DATA in your shell
#export DYSH_DATA="/lma1/teuben/GBT/dysh_data"
export DYSH_DATA="/bigdisk/data/gbt/dysh_data"

dd=("multismallsmall" "multismallbig" "multihugesmall" "multibighuge")
nfile=(8 8 5 3)
#####################
# GBTFITSLOAD
#####################

# skip gbtfitsload because i already ran this bench
for i in $(seq 0 3)
do
    out="benchtest$i.tab"
    #opr="${dd[$i]}.profile.time"
    #oprs="${dd[$i]}.profile.skipflags.time"
    opr="${dd[$i]}.profile"
    oprs="${dd[$i]}.profile.skipflags"
    oprsi="${dd[$i]}.profile.skipflags.useindex"
    #for j in $(seq 1 ${nfile[$i]})
    j=${nfile[$i]}
    #    ../bench_gbtfitsload.py -l 4 -n $j    -k ${dd[$i]} --statslines 50 -m -p  > $opr
        ../bench_gbtfitsload.py -l 4 -n $j --indexthreshold=-1 -s -k ${dd[$i]} --statslines 50 -m -p  > $oprs
        ../bench_gbtfitsload.py -l 4 -n $j --indexthreshold=0 -s -k ${dd[$i]} --statslines 50 -m -p  > $oprsi
        echo "done ${dd[$i]}"
done

#####################
# GETPS
#####################
echo "doing GETPS bench..."
for i in $(seq 0 3)
do
    ../bench_getps.py -s --indexthreshold=-1 -m -k ${dd[$i]} -b fitsio > "getps_bench.fitsio.${dd[$i]}.profile"
    ../bench_getps.py -s --indexthreshold=-1 -m -k ${dd[$i]} -b astropy > "getps_bench.astropy.${dd[$i]}.profile"
    ../bench_getps.py -s --indexthreshold=-1 -m -k ${dd[$i]} > "getps_bench.None.${dd[$i]}.profile"
done
