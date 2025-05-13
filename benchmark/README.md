# Notes for Q8-Q9 benchmarking Dysh
 
The purpose of this  benchmarking is to:

## Q8
- Identify bottlenecks in memory, I/O, or CPU, especially as compared to GBTIDL performance.  
- If possible, indicate (but not implement) potential solutions

## Q9
- Identify and implement solutions to bottlenecks found in Q8.

We wrote a dysh.util.timers.DTime() class that in its simplest form tells us CPU and MEM usage via naming tags:

```
      dt = DTime()
      do_work1()
      dt.tag("work1")
      do)=_work2()
      dt.tag("work2")
      dt.close()
      dt.report()

```
for example here is the output of `bench_getps.py -s -d -t`

```
  name     time VmSize VmRSS  skipflags
             ms  Mbyte  Mbyte           
---------- ----- ------ ------ ---------
      load 101.2 2251.1  311.9      True
   getps1s 193.0 2263.3  326.4      True
   getps1t 614.7 2325.2  390.0      True
   getps2s 298.3 2300.5  365.4      True
   getps2t 137.7 2300.5  365.4      True
   getps3s 185.9 2300.5  365.4      True
   getps3t 132.5 2300.5  365.4      True
   getps4s 184.6 2300.5  365.4      True
   getps4t 132.6 2300.5  365.4      True
```
## 10
- Fit benchmarking data to find pain points

We wrote fitbench.py (still under devel) that allows you to do a linear fit plus offset of one column of the benchmark data (default "time") to up to 3 other columns of the data. e.g.

```
./fitbench.py --files bench*.tab -p nchan nrow nflag
```

It returns the popt, pcov, and np.diag(pcov) from scipy.optimize.curve_fit

# Overall findings

0. Benchmarking is tricky, you measure CPU and MEM usage that does not always make sense. For example we have a case where repeated
   calls to getps() showed unreasonable variations.

1. Overhead in GBTFITSLoad(skipflags=False) - the default - can be very large, notably for ARGUS examples. 9sec vs. 9min were seen.  Our implementation of processing Flag files via Selection object, while convenient, is expensive.  We need to drill down and find out why this is so.

2. Overhead of working on a few scans from a big file, vs. a file containing only those scans seems to suggest that all scans were "used".

3. more to come.


### Notes

  - Not all operations are one-to-one with GBTIDL. For instance GBTIDL cannot calibrate multiple scans at once, whereas dysh can.  We will indicate in tabulated results if a dysh operation has no GBTIDL analog
  - dysh always creates the analog of GBTIDL's index file, so GBTIDL comparisons should be run with no index file.
 - Timing is 'wall clock time', i.e., it includes and kernel/sleep operations.  In python, we are using  *time.perfcounter_ns()*

## Avoiding File caching

On linux files are caching so repeatedly running a benchmark on the same file will improve performance the 2nd time.  
Therefore you have to turn off file caching **each time** you run your benchmark with:

    sudo echo 1 > /proc/sys/vm/drop_caches

if this gives permission denied, open up a root shell, and issue the command in that shell as root:

    sudo su
    sync;sync;sync
    echo 1 > /proc/sys/vm/drop_caches

## disk I/O

Although `hdparm -t` will report a typical I/O speed, in real life this is never achieved. Blocksize of reading affects timing. For dysh we mostly care about read time.

## OMP_NUM_THREADS
The OMP_NUM_THREADS environment variable sets the number of threads to use for  parallel processing.  Peter has in other benchmarks found that changing this from unset to 1 can affect performance.   So you should try your benchmark with both states, e.g. (csh):

    unsetenv OMP_NUM_THREADS

and
    setenv OMP_NUM_THREADS 1

## SDFITS Files

The standard set of SDFITS files to run the benchmark on are:

1. Standard positionswitch example from the notebooks - up to 4 PS on/off scans. Single fits file, 45MB - AGBT05B_047_01
2. L-band edge data in On/Off/On mode - 8.5GB, 70 scans.   NGC2808 is extracted from this, using 9 "Track" scans in on/off/on mode. - AGBT15B_287_19
3. ARGUS edge data in OTF mode - 1.3 GB, NGC0001 in AGBT21B_024_01 (or TBD if NGC5954 (2 strong sources) should be used - AGBT21B_024_20)


### DYSH_DATA

To make it portable accross machines, data not in the $DYSH/testdata directory can be easily used via the `dysh.util.files.dysh_data()` function. Either placed
in $DYSH_DATA/sdfits, under $DYSH_DATA/example-data or $DYSH_DATA/acceptance_testing.   TBD


### Examples

The first example already showed very bizarre behavior when we repeated various forms of `getps()` in a repeated loop. You would expect
times to be the same, but it all depended on if we time-averaged or not.  Times are in ms, on the lma machine at UMD.

```
$ OMP_NUM_THREADS=1 /usr/bin/time bench_getps.py -d -l 10 -t
load     97.532862   Loading this was almost 100ms
getps1s 190.831522   just a getps() to the PSScan is returned
getps1t 672.748917   now .timeaverage() is added to return a Spectrum
getps2s 401.706166
getps2t 137.061235
getps3s 186.651159
getps3t 135.26012
getps4s 186.649279
getps4t 151.799951
getps5s 191.804484
getps5t 136.399024
getps6s 186.610638
getps6t 134.886799
getps7s 184.546092
getps7t 135.14739
getps8s 184.658502
getps8t 133.752485
getps9s 183.76686
getps9t 134.601758
end 0.00999
```
if the -t was not added, only repeated PSScan's were obtained, and all the times was very compatible and about 185ms, especially
the ``getps2s`` stands out at over 400ms.

## NEMO

In NEMO two programs exist that go down to the C level, including an option to use OpenMP, to benchmark a typical 4-phase calibration cycle
of a PS observation (and perhaps more). The work is very similar to bench_sdmath.py, some care should be taken not to make nscan*nchan too small,
or else CPU cache will be used and code will look too fast.  A good compromise seems to be:

      /usr/bin/time sdmath nscan=1000 nchan=100000 iter=10

which typically takes 3.5 on lma, depending on the CPU (e.g. 11.5s on fourier)

### getps

The example='test' in dysh are 8 ON/OFF scans, 352 rows, thus nscan=88 for nchan=32768.


      python bench_getps.py  -d -s -l 10 
      -> 600 ms

      bench_sdmath.py -d -s 88 -n 32768 -m -1 -l 10
      -> 18 ms

      /usr/bin/time sdmath 88 32768 1000
      -> 10.6 ms

Also perhaps to note is that bench_sdmath runs faster in float32/float64 (60-90% in two cases), whereas the C program
only ran 20% faster. Odd.

        # getps, timing in ms
        # with 4 scans     with 1 scan
 
        -n      270          100
        .       600          180
        -t      950          315

### gbtfitsload

## TODO

- API args= to **kwargs ?
- always do CPU and MEM, currently MEM is done via args.memory=True
- the command
       python -m cProfile  ./bench_getps.py -s -d -t -l 1 |less
  gives different results from 
       ./bench_getps.py -s -d -t -l 1 -p
  Most notably, it claims in the former that loading gb20mfitsload.py took 2.7sec
