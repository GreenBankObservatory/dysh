#!/bin/bash

# Script to run various benchmarks
OMP_NUM_THREADS=1

export DYSH_DATA="/lma1/teuben/GBT/dysh_data"

dd=("multismallsmall" "multismallbig" "multihugesmall" "multibighuge")
nfile=(8 8 5 3)
#####################
# GBTFITSLOAD
#####################

for i in $(seq 0 3)
do
    out="benchtest$i.tab"
    opr="${dd[$i]}.profile"
    oprs="${dd[$i]}.profile.skipflags"
    #for j in $(seq 1 ${nfile[$i]})
    j=${nfile[$i]}
    #do
        ./bench_gbtfitsload.py -d -l 4 -n $j  -k ${dd[$i]} --statslines 50 -m -p  > $opr
        ./bench_gbtfitsload.py -d -l 4 -n $j -s -k ${dd[$i]}  --statslines 50 -m  -p  > $oprs
        echo "done ${dd[$i]}"
    #done
done

exit
