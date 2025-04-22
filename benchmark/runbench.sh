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
    opr="$out.pr"
    oprs="$out.prs"
    ./bench_gbtfitsload.py -d -l 4 -n ${nfile[$i]} -s -k ${dd[$i]} -a -o $out -p -m  > $opr
    ./bench_gbtfitsload.py -d -l 4 -n ${nfile[$i]} -k ${dd[$i]} -a -o $out -p -m >  $oprs
done

exit

datakey="multismallsmall"
OUTTAB="benchtest1.tab"
for num in $(seq 1 8);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
datakey="multismallbig"
OUTTAB="benchtest2.tab"
for num in $(seq 1 8);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
datakey="multihugesmall"
OUTTAB="benchtest3.tab"
for num in $(seq 1 5);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  

datakey="multibighuge"
OUTTAB="benchtest4.tab"
for num in $(seq 1 3);
do
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -s -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.prs"
   ./bench_gbtfitsload.py -d -l 4 -n ${num}  -k ${datakey} -a -o ${OUTTAB} -p -m > "${OUTTAB}.pr"
done
# now just print the final table
./bench_gbtfitsload.py -j -o ${OUTTAB} 
  
