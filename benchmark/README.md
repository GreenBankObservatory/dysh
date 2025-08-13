# Notes for Q8-Q9 benchmarking Dysh
 
The purpose of this  benchmarking is to:

## Q8
- Identify bottlenecks in memory, I/O, or CPU, especially as compared to GBTIDL performance.  
- If possible, indicate (but not implement) potential solutions

THe benchmarks created are:


* bench_datawrite.py: Test writing of data performance

* bench_edge.py: Test performance of data writing using the EDGE data.  This is not a standard benchmark

* bench_gbtfitsload.py: Test performance of loading SDFITS files with and without flagging

* bench_getps.py:  Test performance of getps, calibration, and averaging.

* bench_otf.py:  Test performance of OTF calibration

* bench_sdmath.py: Test performance of numpy equivalent of getps.  Can be used as a baseline of sorts of the fastest we could do these operations.

* bench_skel.py:  Skeleton (template) from which new benchmarks can be created.

## Q9
- Identify and implement solutions to bottlenecks found in Q8.


# Methodology

We wrote a `dysh.util.timers.DTime()` class that in its simplest form tells us CPU and MEM usage via naming tags:

```
      dt = DTime()
      do_work1()
      dt.tag("work1")
      do_work2()
      dt.tag("work2")
      dt.close()
      dt.report()

```
for example here is the output of `./bench_getps.py -s -d -t` on our favorite `lma` machine:

```
  name     time VmSize VmRSS  skipflags
             ms  Mbyte  Mbyte           
---------- ----- ------ ------ ---------
   load 398.1 1684.7 315.5      True
getps1t 763.2 1795.8 438.5      True
getps2t 414.7 1790.9 433.6      True
getps3t 510.4 1764.9 408.0      True
getps4t 439.4 1764.9 408.0      True
 report   0.0 1764.9 408.0      True
final 2.525809781  sec

```

- Fit benchmarking data to find pain points

We wrote fitbench.py (still under devel) that allows you to do a linear fit plus offset of one column of
the benchmark data (default "time") to up to 3 other columns of the data. e.g.

```
./fitbench.py --files bench*.tab -p nchan nrow nflag
```

It returns the popt, pcov, and np.diag(pcov) from scipy.optimize.curve_fit

# Overall findings

1. Benchmarking is tricky, you measure CPU and MEM usage that does not
   always make sense. For example we have a case where repeated calls
   to getps() showed unreasonable variations.
   
2. Overhead in `GBTFITSLoad(skipflags=False)` - the default - can be
   very large, notably for ARGUS examples. 9sec vs. 9min were seen.
   Our implementation of processing Flag files via Selection object,
   while convenient, is expensive.

3. Overhead of working on a few scans from a big file, vs. a file
   containing only those scans seems to suggest that all scans were
   "used".  (#678 memory leak)

4. There is sometimes an extra overhead on the first of many loops, as
   the `getps` shows above. It is unclear where this "setup" times comes
   from possibly.  Possibly python learning how to set up a class. But
   not something we can work around. Same with I/O, on the second
   benchmark run, the I/O overhead in the `load` step can be significantly
   smaller. 


More details are in `q8stats/README.md`.

### Notes

  - Not all operations are one-to-one with GBTIDL. For instance, GBTIDL cannot calibrate multiple scans at once, whereas dysh can.  
  - dysh always creates the analog of GBTIDL's index file, so GBTIDL comparisons should be run with no index file.
  - Timing is 'wall clock time', i.e., it includes and kernel/sleep operations.  In python, we are using  *time.perfcounter_ns()*

## Avoiding File caching

On linux files are caching so repeatedly running a benchmark on the same file will improve performance the 2nd time.  
Therefore it is best if you turn off file caching **each time** you run your benchmark with:

    sudo echo 1 > /proc/sys/vm/drop_caches

if this gives permission denied, open up a root shell, and issue the command in that shell as root:

    sudo su
    sync;sync;sync
    echo 1 > /proc/sys/vm/drop_caches
or
    echo 1 | sudo tee /proc/sys/vm/drop_caches

Example  reading a 8GB file on a laptop:

    ./bench_datawrite.py --source 1 --writedata sb
    -------- ------- ------- ------ ------ --------- ------- ----- ----- --- --- ---- ------ --------- ------
    init     4.8  1627.2  276.4      0       0.0     0.0     0     0   0   0    0      0     False      0
    load 44366.5 11879.2 6138.2      2   3963.61 7927.22 32768 60192   8   1    2     28     False     -1
    getps 2337.1 12042.6 6341.4      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
    -------- ------ ------- ------ ------ --------- ------- ----- ----- --- --- ---- ------ --------- ------
    init     2.1  1627.2  277.7      0       0.0     0.0     0     0   0   0    0      0     False      0
    load  2752.3 11879.2 6140.8      2   3963.61 7927.22 32768 60192   8   1    2     28     False     -1
    getps 2350.3 12042.6 6342.1      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1


## disk I/O

Although `hdparm -t` will report a typical I/O speed, in real life this is never achieved.
Blocksize of reading affects timing. For dysh we mostly care about read time.

## OMP_NUM_THREADS

The OMP_NUM_THREADS environment variable sets the number of threads to
use for parallel processing.  Peter has in other benchmarks found that
changing this from unset to 1 can affect performance.  So you should
try your benchmark with both states, e.g. (csh):

    unsetenv OMP_NUM_THREADS

and
    setenv OMP_NUM_THREADS 1

In the future this may be resolved when the GIL (Global Interpreter Lock) has been
resolved.

## SDFITS Files

The standard set of SDFITS files to run the benchmark on are:

1. Standard L-band positionswitch example from the notebooks - up to 4 PS on/off scans. Single fits file, 45MB - AGBT05B_047_01.
   This is **example=*getps"**

2. L-band edge data in On/Off/On mode - 8.5GB, 70 scans.   NGC2808 is extracted from this, using 9 "Track" scans in on/off/on mode. - AGBT15B_287_19.
   This should be **example=*getps3"**   ???

3. ARGUS edge data in OTF mode - 1.3 GB, NGC0001 in AGBT21B_024_01 (or TBD if NGC5954 (2 strong sources) should be used - AGBT21B_024_20).
   This is **example=*otf1"**


### DYSH_DATA

To make it portable accross machines, data not in the $DYSH/testdata
directory can be easily used via the `dysh.util.files.dysh_data()`
function. Either placed in $DYSH_DATA/sdfits, under
$DYSH_DATA/example-data or $DYSH_DATA/acceptance_testing.


### Summary of tests and Results

The full report is in `q8stats/README.md`.


#### Example: Strange behavior when time-averaging

The first example already showed very bizarre behavior when we repeated various forms of `getps()` in a repeated loop. You would expect
times to be the same, but it all depended on if we time-averaged or not.  Times are in ms, on the lma machine at UMD.  This may be related to issue #583.

```
$ OMP_NUM_THREADS=1 /usr/bin/time bench_getps.py -d -l 10 -t
...
  name    time  VmSize VmRSS skipflags
           ms   Mbyte  Mbyte          
-------- ------ ------ ----- ---------
    load  220.0  964.8 325.9     False
 getps1t 1163.6 1066.3 439.2     False
 getps2t  673.6 1063.9 436.8     False
 getps3t  683.1 1064.6 437.5     False
 getps4t  671.7 1064.6 437.5     False
 getps5t  674.5 1065.5 438.5     False
 getps6t  674.5 1065.2 438.0     False
 getps7t  681.0 1064.2 437.0     False
 getps8t  675.9 1064.6 437.5     False
 getps9t  674.0 1064.6 437.5     False
getps10t  675.8 1064.6 437.5     False
  report    0.1 1064.6 437.5     False

```
if the -t was not added, the first measurment did not stand out as much. x

## Disk I/O

There are three types of disk I/O:  the sdfits, the scanblock and the spectrum.

Few comments here:

1. Disk caching clearly plays a role, especially on machines with modest memory.
2. Python 


```
# writing sp and sb have similar properties


  name    time   VmSize VmRSS  #files file_size totsize nchan  nrow nIF nFd nPol #flags skipflags nwrite
           ms    Mbyte  Mbyte           Mbyte    Mbyte                                                  
------- ------- ------- ------ ------ --------- ------- ----- ----- --- --- ---- ------ --------- ------
   load 43578.5 11879.2 6137.1      2   3963.61 7927.22 32768 60192   8   1    2     28     False     -1
  getps  2402.2 12047.3 6345.0      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite   483.3 12126.5 6424.4      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite   275.8 12070.5 6369.1      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite   154.4 12070.5 6369.1      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1

------- ------ ------- ------ ------ --------- ------- ----- ----- --- --- ---- ------ --------- ------
   load 2754.4 11878.2 6140.7      2   3963.61 7927.22 32768 60192   8   1    2     28     False     -1
  getps 2331.0 12042.6 6342.6      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite  444.2 12126.5 6426.8      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite  282.4 12070.5 6371.6      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
spwrite  156.8 12070.5 6371.6      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1

--------- ------- ------- ------- ------ --------- ------- ----- ----- --- --- ---- ------ --------- ------
     load  2856.7 11879.2  6248.8      2   3963.61 7927.22 32768 60192   8   1    2     28     False     -1
    getps  2404.9 12047.3  6455.5      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
sdfwrite1  6871.7 26296.3 24184.7      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1
sdfwrite2 79575.9 40552.5 16230.8      2   3963.61 7927.22 32768 60192   8   1    2     28     False      1

```

       a.   First time `load` was about 43sec (the SDFITS file is 7GB), after that more like 3 sec.
            This is bad for science, as only there the first one counts.
       b.   The `spwrite` stage seems to also suffer from a speedup as it was re-executed for the benchmark.
            This is again bad for science, as only there the first one counts. In fact, only on the 3rd
            execution did it flatten out as about 155 ms.
       c.   As already seen in the `getps` benchark, there was also a modest re-execution penalty for `getps`,
            but quite noticable when time-averaging was added. (80%)
       d.   During `sdfwrite` a large block of memory was not freed (issue #678) and as a result all subsequent
            runs on a 32GB machine got slowed down.   On a 512 GB machine this was not noticable until 100 loops.


```



For science, we would count the benchmark as 43.6 + 2.4 + 0.5  = 46.5, taking the first-time values,
but if one would take the repeated benchmark values this would count as 2.8 + 2.4 + 0.2 = 5.4 sec. Quite a difference!

## Math and sdmath

The math is relatively simple, and very parallizable:

       Tsys = Tc  <cold> / <hot-cold>       (1)
       Ta   = Tsys  (on-off)/off            (2)

First of all, this implies a Tsys for the ON and one for the OFF.   Which one to use in (2) ?

Secondly, what is the ON and the OFF really? Can one use the calon + caloff average as in:

       (on_calon+on_caloff) / (off_calon+off_caloff) -1


### NEMO

In NEMO two programs exist that go down to the C level, including an option to use OpenMP, to benchmark a typical 4-phase calibration cycle
of a PS observation (and perhaps more). The work is very similar to bench_sdmath.py, some care should be taken not to make nscan*nchan too small,
or else CPU cache will be used and code will look too fast.  A good compromise seems to be:

      /usr/bin/time sdmath nscan=1000 nchan=100000 iter=10

which typically takes 3.5 on lma, depending on the CPU (e.g. 11.5s on fourier)

## TODO

- API args= to **kwargs ?
- always do CPU and MEM, currently MEM is done via args.memory=True
- the command
       python -m cProfile  ./bench_getps.py -s -d -t -l 1 |less
  gives different results from 
       ./bench_getps.py -s -d -t -l 1 -p
  Most notably, it claims in the former that loading gb20mfitsload.py took 2.7sec
