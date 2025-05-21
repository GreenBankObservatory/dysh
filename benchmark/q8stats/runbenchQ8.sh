#!/bin/bash

# Script to run various benchmarks
OMP_NUM_THREADS=1

export DYSH_DATA="/lma1/teuben/GBT/dysh_data"

dd=("multismallsmall" "multismallbig" "multihugesmall" "multibighuge")
nfile=(8 8 5 3)
#####################
# GBTFITSLOAD
#####################

# skip gbtfitsload because i already ran this bench
for i in $(seq 0 3)
do
    out="benchtest$i.tab"
    opr="${dd[$i]}.profile.time"
    oprs="${dd[$i]}.profile.skipflags.time"
    #for j in $(seq 1 ${nfile[$i]})
    j=${nfile[$i]}
    #do
        ../bench_gbtfitsload.py -d -l 4 -n $j  -k ${dd[$i]} --statslines 50 -m -p  -x time > $opr
        ../bench_gbtfitsload.py -d -l 4 -n $j -s -k ${dd[$i]}  --statslines 50 -m  -p  -x time > $oprs
        echo "done ${dd[$i]}"
    #done
done

#####################
# GETPS
#####################
echo "doing GETPS bench..."
../bench_getps.py -d -s -t --statslines 50 -m -p -x time > getps_bench_ta.profile.time
../bench_getps.py -d -s -t --statslines 50 -m -p > getps_bench_ta.profile
../bench_getps.py -d -s --statslines 50 -m -p -x time > getps_bench.profile.time
../bench_getps.py -d -s --statslines 50 -m -p > getps_bench.profile

#####################
# SCAN WRITE
#####################
echo "doing Scanblock write bench..."
../bench_datawrite.py --writedata sb -l 10 --statslines 25 -m -p > sbwrite_bench.profile
exit
