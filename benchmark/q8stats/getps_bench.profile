using Namespace(dobench=True, key='test1', timeaverage=False, nocalibrate=False, loop=4, skipflags=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='cumulative', memory=True, quit=False)
Loading  /lma1/teuben/GBT/dysh_data/example_data/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits
STATS: {'nrows': 352, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 2, 'intnum': 11, 'sig': 1, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: getps
  name   time  VmSize VmRSS skipflags
          ms   Mbyte  Mbyte          
------- ------ ------ ----- ---------
   load  243.1 3156.0 247.3      True
getps1s 1336.7 3198.4 301.4      True
getps2s 1322.9 3229.0 331.4      True
getps3s 1326.5 3229.7 332.1      True
getps4s 1324.5 3240.7 343.1      True
 report    0.1 3240.7 343.1      True
         8488540 function calls (8153134 primitive calls) in 5.560 seconds

   Ordered by: cumulative time
   List reduced from 1871 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.001    0.000    5.309    1.327 /lma1/mpound/dysh/src/dysh/log.py:277(wrapper)
        4    0.002    0.001    5.306    1.326 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:1339(getps)
       16    0.001    0.000    4.764    0.298 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1157(__init__)
       16    0.000    0.000    2.705    0.169 /lma1/mpound/dysh/src/dysh/spectra/scan.py:172(_finish_initialization)
      352    0.052    0.000    2.174    0.006 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1338(exposure)
      233    0.010    0.000    2.028    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      166    0.004    0.000    1.921    0.012 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
     1440    0.004    0.000    1.704    0.001 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:124(index)
    11620    0.185    0.000    1.700    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
     4651    0.019    0.000    1.656    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4062(__getitem__)
       64    0.001    0.000    1.583    0.025 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:380(rawspectra)
     1537    0.006    0.000    1.489    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4130(_getitem_bool_array)
25077/17358    0.019    0.000    1.460    0.000 {method 'view' of 'numpy.ndarray' objects}
       16    0.007    0.000    1.416    0.088 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1276(calibrate)
       83    0.005    0.000    1.371    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
     1480    0.004    0.000    1.349    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6662(copy)
     1517    0.008    0.000    1.329    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:557(copy)
       83    0.034    0.000    1.273    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
       16    0.001    0.000    1.084    0.068 /lma1/mpound/dysh/src/dysh/spectra/scan.py:489(_add_calibration_meta)
     1458    0.082    0.000    0.911    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:1782(_consolidate_inplace)
378610/204307    0.110    0.000    0.860    0.000 {built-in method builtins.setattr}
   185920    0.177    0.000    0.768    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
     1446    0.031    0.000    0.720    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:2259(_consolidate)
       82    0.004    0.000    0.568    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1512(_init_from_coldefs)
     5740    0.006    0.000    0.557    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1638(_copy_column)
     5740    0.010    0.000    0.548    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:960(copy)
     2900    0.010    0.000    0.534    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexing.py:1176(__getitem__)
        4    0.001    0.000    0.495    0.124 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:977(_common_selection)
    11620    0.017    0.000    0.485    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:871(name)
       18    0.000    0.000    0.449    0.025 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:544(nchan)
       18    0.000    0.000    0.449    0.025 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:404(rawspectrum)
     2880    0.009    0.000    0.447    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexing.py:1719(_getitem_axis)
     1552    0.005    0.000    0.410    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:4142(_take_with_is_copy)
     1561    0.008    0.000    0.401    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
     1565    0.085    0.000    0.397    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:317(apply)
     1561    0.005    0.000    0.368    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:869(take)
     1440    0.001    0.000    0.363    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexing.py:1696(_get_list_axis)
    10120    0.157    0.000    0.359    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:2276(_merge_blocks)
     6131    0.011    0.000    0.313    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:161(__init__)
     2022    0.018    0.000    0.305    0.000 {built-in method builtins.sorted}
   141724    0.042    0.000    0.300    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:2264(<lambda>)
    25475    0.016    0.000    0.294    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/config/configuration.py:333(__get__)
     1561    0.020    0.000    0.293    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:623(reindex_indexer)
    25475    0.069    0.000    0.278    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/config/configuration.py:442(__call__)
    70862    0.081    0.000    0.258    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/blocks.py:225(_consolidate_key)
58256/254    0.055    0.000    0.246    0.001 /lma1/mpound/Anaconda3-2024.10-1/lib/python3.12/copy.py:118(deepcopy)
    12223    0.020    0.000    0.246    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/blocks.py:1287(take_nd)
        4    0.000    0.000    0.241    0.060 /lma1/mpound/dysh/src/dysh/util/selection.py:871(__deepcopy__)
1097863/1094378    0.189    0.000    0.231    0.000 {built-in method builtins.isinstance}
        4    0.000    0.000    0.223    0.056 /lma1/mpound/dysh/src/dysh/util/selection.py:838(_select_from_mixed_kwargs)


final 5.5537966370000005  sec
