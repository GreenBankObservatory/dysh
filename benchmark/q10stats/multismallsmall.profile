Loading not more than 8 from /lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT20B_336_01/AGBT20B_336_01.raw.vegas
Will load 8 of 8 files. FITS size per file 56.18MB, Flag lines 7
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 8 FITS files
Found ostype= linux
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 8 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 8 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 8 FITS files
# Dysh Benchmark: GBTFITSLoad
 name  time  VmSize VmRSS  #files file_size totsize nchan  nrow nIF nFd nPol #flags skipflags
        ms   Mbyte  Mbyte           Mbyte    Mbyte                                           
----- ------ ------ ------ ------ --------- ------- ----- ----- --- --- ---- ------ ---------
load1 7500.7 1653.1 1242.6      8     56.18  449.44  1024 95136   4   2    2      7     False
load2 7532.5 2101.7 1691.1      8     56.18  449.44  1024 95136   4   2    2      7     False
load3 7480.0 2193.8 1783.2      8     56.18  449.44  1024 95136   4   2    2      7     False
load4 7453.9 2195.3 1784.4      8     56.18  449.44  1024 95136   4   2    2      7     False
         51375072 function calls (50914623 primitive calls) in 29.783 seconds

   Ordered by: cumulative time
   List reduced from 1819 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   29.088    7.272 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.001    0.000   29.086    7.271 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
       32    0.001    0.000   20.405    0.638 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
       32    0.033    0.001   19.341    0.604 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
     1508    2.609    0.002   10.235    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       32    0.000    0.000    9.682    0.303 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       32    0.025    0.001    9.680    0.303 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       32    0.000    0.000    9.655    0.302 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       32    0.003    0.000    9.577    0.299 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      736    0.002    0.000    9.500    0.013 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
        4    0.231    0.058    8.506    2.127 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:892(_create_index_if_needed)
       32    0.000    0.000    6.907    0.216 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:505(drop_fields)
       32    2.266    0.071    6.904    0.216 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:32(recursive_fill_fields)
      772    0.024    0.000    6.465    0.008 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
     2408    0.008    0.000    5.615    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      736    0.005    0.000    4.788    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
     2372    0.012    0.000    4.632    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:692(field)
     2372    0.010    0.000    4.544    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:883(_convert_other)
      736    0.002    0.000    4.502    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/util.py:254(decode_ascii)
      736    0.053    0.000    4.498    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:509(decode)
      736    0.002    0.000    4.476    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
      736    0.001    0.000    4.391    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
  8752512    1.702    0.000    3.128    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
      736    3.128    0.004    3.128    0.004 {built-in method numpy._core._multiarray_umath._vec_string}
       32    0.003    0.000    2.900    0.091 /lma1/mpound/dysh/src/dysh/util/selection.py:1281(read)
  8752512    1.661    0.000    2.832    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
        4    0.068    0.017    2.591    0.648 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:970(_construct_integration_number)
    11904    0.041    0.000    2.373    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
    11904    0.012    0.000    2.324    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1149(__iter__)
      136    0.004    0.000    2.250    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
       32    0.001    0.000    2.153    0.067 /lma1/mpound/dysh/src/dysh/util/selection.py:1124(flag)
        4    0.011    0.003    2.081    0.520 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:944(_construct_procedure)
       36    0.001    0.000    1.965    0.055 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:834(split)
    11896    0.678    0.000    1.899    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1180(_chop)
       32    0.001    0.000    1.842    0.058 /lma1/mpound/dysh/src/dysh/util/selection.py:534(_base_select)
    10064    0.166    0.000    1.637    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
47396/23364    0.054    0.000    1.548    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/decorators.py:842(__get__)
       32    0.000    0.000    1.480    0.046 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:404(data)
       32    0.000    0.000    1.459    0.046 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:174(_get_tbdata)
  8752527    1.426    0.000    1.426    0.000 {method 'decode' of 'bytes' objects}
55714/49118    0.032    0.000    1.388    0.000 {method 'view' of 'numpy.ndarray' objects}
      736    0.015    0.000    1.316    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:98(_to_bytes_or_str_array)
469081/318116    0.142    0.000    1.296    0.000 {built-in method builtins.setattr}
28288/26712    1.215    0.000    1.217    0.000 {built-in method numpy.asarray}
       68    0.004    0.000    1.214    0.018 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
  8815192    1.186    0.000    1.186    0.000 {method 'strip' of 'str' objects}
       96    0.014    0.000    1.137    0.012 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3575(sort)
       36    0.000    0.000    1.136    0.032 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:327(_str_split)
       68    0.030    0.000    1.126    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
    37667    0.115    0.000    1.036    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:704(__array_finalize__)
