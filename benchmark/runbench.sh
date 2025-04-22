#!/bin/bash

# Script to run various benchmarks
OMP_NUM_THREADS=1

export DYSH_DATA="/lma1/teuben/GBT/dysh_data"

#####################
# GBTFITSLOAD
#####################

datakey="multismallsmall"
OUTTAB="benchtest1.tab"
for num in $(seq 1 8);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
datakey="multismallbig"
OUTTAB="benchtest2.tab"
for num in $(seq 1 8);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
datakey="multihugesmall"
OUTTAB="benchtest3.tab"
for num in $(seq 1 5);
do
   #echo  ./bench_gbtfitsload.py -d -l 4 -n ${num} -a -o ${OUTTAB}  -s -k ${datakey}
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  

datakey="multibighuge"
OUTTAB="benchtest4.tab"
for num in $(seq 1 3);
do
   #echo  ./bench_gbtfitsload.py -d -l 4 -n ${num} -a -o ${OUTTAB}  -s -k ${datakey}
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
