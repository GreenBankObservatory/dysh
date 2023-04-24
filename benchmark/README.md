# Notes for benchmarking Dysh

The purpose of  benchmarking is to test performance  of 

- opening/loading an SDFITS file/hdu 
- optionally creating an index of some sort (GBTIDL or pandas)
- creating spectra from each DATA row of the SDFITS table.
- remove baselines of order 1, 2, and 3 from each spectrum.
 
Timing is 'wall clock time', i.e., it includes and kernel/sleep operations.  In python, we are using  *time.perfcounter_ns()*

## Avoiding File caching
On linux files are caching so repeatedly running a benchmark on the same file will improve performance the 2nd time.  Therefore you have to turn off file caching **each time** you run your benchmark with:

    sudo echo 1 > /proc/sys/vm/drop_caches

## OMP_NUM_THREADS
The OMP_NUM_THREADS environment variable sets the number of threads to use for  parallel processing.  Peter has in other benchmarks found that changing this from unset to 1 can affect performance.   So you should try your benchmark with both states, e.g. (csh):

    unsetenv OMP_NUM_THREADS

and
    setenv OMP_NUM_THREADS 1

## SDFITS Files 

The standard set of SDFITS files to run the benchmark on are in /lma1/mpound/GBT/examples/ on LMA machines and /home/scratch/mpound/examples in GBO machines.  The files are (in ascending order of size):

- rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.E.fits
- misc/IC1481.fits
- misc/W3OH.fits
- misc/ngc5291.fits
- onoff-L/data/TGBT21A_501_11.raw.vegas.fits
- mapping-L/data/TGBT17A_506_11.raw.vegas/
TGBT17A_506_11.raw.vegas.A.fits
- nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.F.fits 
- mixed-fs-ps/data/AGBT16B_225_05/AGBT16B_225_05.raw.vegas/AGBT16B_225_05.raw.vegas.B.fits

## Code details
1. Benchmark code is in the subdirectory *benchmark*
2. The python benchmark code is *revisedstructure.py* and is invoked with *benchmark_py.csh*
    -  To run it, you have to have dysh installed.  This is done with hatch and I haven't written up the notes for that yet!
    - It has various options to turn on/off pieces of the benchmark, see *revisedstructure.py -h*
 
## Output
Output should be into a (IPAC-formatted) table with the following columns:

- File - file name
- Size - file size MB
- N_hdu - number of HDU in the file
- HDU - number of the HDU for which this row contains the benchmark
- N_rows - number of rows in this HDU
- N_chan - number of chan per spectum in this HDU  
- Load - time to load/read the SDFITs file, in ms .  e.g. astropy.io.fits.open()
- Index - time to create index (e.g. pandas, GBTIDL) in ms
- Create_Obslocks - time to create Spectrum object for every row of HDU, in ms
- Baseline_N for N=1,2,3  - time to remove baseline of order N from every row of HDU, in ms
- Total - total time taken in ms (sum of other columns)

# Plotting
The python script *makeplots.py* will make some standard  plots given the IPAC formatted table.   *makeplots.py -h* to see the options. e.g.

    ./makeplots.py -f tab_omp=1.out -t "Timing for Load/Obsblocks/Baseline (OMP=1)" -o t1.png
    ./makeplots.py -f tab_omp_unset.out -t "Timing for Load/Obsblocks/Baseline (OMP unset)" -o t2.png

