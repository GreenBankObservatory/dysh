#

t1=$(./bench_getps.py -d  -l 9                            | grep ^getps | awk '{sum += $2} END{print sum/NR}')
t2=$(./bench_sdmath.py -d --nscan 4400 --nchan 32768 -l 9 | grep ^math2 | awk '{sum += $2} END{print sum/NR}')

echo $t1 $t2 | awk '{print $1,$2/$1}'
