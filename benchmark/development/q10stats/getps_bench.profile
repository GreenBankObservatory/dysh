using Namespace(dobench=True, key='test1', timeaverage=False, nocalibrate=False, loop=4, skipflags=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='cumulative', memory=True, quit=False)
AGBT05B_047_01.raw.acs.fits already downloaded
Loading  AGBT05B_047_01.raw.acs.fits
STATS: {'nrows': 352, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 2, 'intnum': 11, 'sig': 1, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: getps
  name   time VmSize VmRSS skipflags
          ms  Mbyte  Mbyte          
------- ----- ------ ----- ---------
   load 240.5 3156.1 249.1      True
getps1s 833.5 3198.5 303.0      True
getps2s 839.6 3228.0 333.0      True
getps3s 838.3 3228.0 333.0      True
getps4s 836.6 3238.7 343.5      True
 report   0.1 3238.7 343.5      True
         5958385 function calls (5647472 primitive calls) in 3.594 seconds

   Ordered by: cumulative time
   List reduced from 1943 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.002    0.000    3.345    0.836 /lma1/mpound/dysh/src/dysh/log.py:282(wrapper)
        4    0.003    0.001    3.340    0.835 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:1388(getps)
       16    0.001    0.000    2.752    0.172 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1499(__init__)
      233    0.011    0.000    2.002    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      166    0.004    0.000    1.894    0.011 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
    11620    0.185    0.000    1.676    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
       64    0.001    0.000    1.567    0.024 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:380(rawspectra)
22677/14958    0.016    0.000    1.436    0.000 {method 'view' of 'numpy.ndarray' objects}
       83    0.005    0.000    1.352    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
       83    0.034    0.000    1.255    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
378610/204307    0.105    0.000    0.846    0.000 {built-in method builtins.setattr}
   185920    0.175    0.000    0.756    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
       16    0.001    0.000    0.723    0.045 /lma1/mpound/dysh/src/dysh/spectra/scan.py:301(_finish_initialization)
       82    0.004    0.000    0.560    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1512(_init_from_coldefs)
     5740    0.006    0.000    0.549    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1638(_copy_column)
     5740    0.010    0.000    0.540    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:960(copy)
        4    0.001    0.000    0.509    0.127 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:984(_common_selection)
    11620    0.016    0.000    0.482    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:871(name)
       18    0.000    0.000    0.441    0.024 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:544(nchan)
       18    0.000    0.000    0.440    0.024 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:404(rawspectrum)
     6131    0.011    0.000    0.312    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:161(__init__)
       16    0.007    0.000    0.308    0.019 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1631(calibrate)
    25475    0.017    0.000    0.291    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/config/configuration.py:333(__get__)
    25475    0.068    0.000    0.274    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/config/configuration.py:442(__call__)
      999    0.004    0.000    0.272    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4062(__getitem__)
58248/254    0.054    0.000    0.244    0.001 /lma1/mpound/Anaconda3-2024.10-1/lib/python3.12/copy.py:118(deepcopy)
        4    0.000    0.000    0.238    0.060 /lma1/mpound/dysh/src/dysh/util/selection.py:873(__deepcopy__)
      305    0.001    0.000    0.236    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4130(_getitem_bool_array)
        4    0.000    0.000    0.214    0.054 /lma1/mpound/dysh/src/dysh/util/selection.py:840(_select_from_mixed_kwargs)
        9    0.000    0.000    0.214    0.024 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
    11690    0.023    0.000    0.213    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:996(_verify_keywords)
        4    0.000    0.000    0.212    0.053 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
   174300    0.211    0.000    0.211    0.000 {built-in method builtins.locals}
       16    0.002    0.000    0.211    0.013 /lma1/mpound/dysh/src/dysh/spectra/scan.py:711(_make_meta)
        1    0.000    0.000    0.211    0.211 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
      160    0.001    0.000    0.200    0.001 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:124(index)
     1604    0.006    0.000    0.185    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexing.py:1176(__getitem__)
    11690    0.013    0.000    0.184    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1288(_determine_formats)
      232    0.001    0.000    0.163    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:557(copy)
      196    0.001    0.000    0.162    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6662(copy)
     5880    0.021    0.000    0.155    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:305(value)
35775/35064    0.042    0.000    0.154    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/decorators.py:842(__get__)
   454/30    0.001    0.000    0.153    0.005 /lma1/mpound/Anaconda3-2024.10-1/lib/python3.12/copy.py:247(_reconstruct)
      644    0.011    0.000    0.151    0.000 {method '__deepcopy__' of 'numpy.ndarray' objects}
    12479    0.015    0.000    0.150    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:284(value)
     5880    0.002    0.000    0.147    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:212(__str__)
 1720/104    0.001    0.000    0.147    0.001 /lma1/mpound/Anaconda3-2024.10-1/lib/python3.12/copy.py:252(<genexpr>)
        8    0.000    0.000    0.146    0.018 /lma1/mpound/Anaconda3-2024.10-1/lib/python3.12/copy.py:200(_deepcopy_tuple)
     5880    0.004    0.000    0.145    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/card.py:521(image)
744120/743052    0.121    0.000    0.143    0.000 {built-in method builtins.isinstance}


final 3.588445606  sec
