#!/bin/bash

# Script to run various benchmarks
OMP_NUM_THREADS=1


#####################
# GBTFITSLOAD
#####################

OUTTAB="benchtest.tab"
for num in $(seq 1 8);
do
   echo  ./bench_gbtfitsload.py -d -l 4 -n ${num} -a -o ${OUTTAB} 
   ./bench_gbtfitsload.py -d -l 4 -n ${num} -a -o ${OUTTAB} 
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
#####################


