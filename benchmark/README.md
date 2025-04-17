# Notes for Q8-Q9 benchmarking Dysh
 
The purpose of this  benchmarking is to:

## Q8
- Identify bottlenecks in memory, I/O, or CPU, especially as compared to GBTIDL performance.  
- If possible, indicate (but not implement) potential solutions

## Q9
- Identify and implement solutions to bottlenecks found in Q8.

### Notes

  - Not all operations are one-to-one with GBTIDL. For instance GBTIDL cannot calibrate multiple scans at once, whereas dysh can.  We will indicate in tabulated results if a dysh operation has no GBTIDL analog
  - dysh always creates the analog of GBTIDL's index file, so GBTIDL comparisons should be run with no index file.
 - Timing is 'wall clock time', i.e., it includes and kernel/sleep operations.  In python, we are using  *time.perfcounter_ns()*

## Avoiding File caching
On linux files are caching so repeatedly running a benchmark on the same file will improve performance the 2nd time.  Therefore you have to turn off file caching **each time** you run your benchmark with:

    sudo echo 1 > /proc/sys/vm/drop_caches

## OMP_NUM_THREADS
The OMP_NUM_THREADS environment variable sets the number of threads to use for  parallel processing.  Peter has in other benchmarks found that changing this from unset to 1 can affect performance.   So you should try your benchmark with both states, e.g. (csh):

    unsetenv OMP_NUM_THREADS

and
    setenv OMP_NUM_THREADS 1

## SDFITS Files

The standard set of SDFITS files to run the benchmark on are ...
