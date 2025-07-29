Loading not more than 8 from /lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT23A_432_03/AGBT23A_432_03.raw.vegas
Will load 8 of 8 files. FITS size per file 66.72MB, Flag lines 20
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
 name  time  VmSize VmRSS  #files file_size totsize nchan nrow nIF nFd nPol #flags skipflags
        ms   Mbyte  Mbyte           Mbyte    Mbyte                                          
----- ------ ------ ------ ------ --------- ------- ----- ---- --- --- ---- ------ ---------
load1 9814.3 1302.7  894.9      8     66.72  533.76 16384 8064   4   2    2     20     False
load2 9939.7 1360.7  952.3      8     66.72  533.76 16384 8064   4   2    2     20     False
load3 9836.4 1360.9  953.1      8     66.72  533.76 16384 8064   4   2    2     20     False
load4 9835.4 1477.9 1070.1      8     66.72  533.76 16384 8064   4   2    2     20     False
         64040004 function calls (61814879 primitive calls) in 39.400 seconds

   Ordered by: cumulative time
   List reduced from 1819 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   39.195    9.799 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.002    0.000   39.193    9.798 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
        4    0.014    0.004   34.033    8.508 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:892(_create_index_if_needed)
       32    0.037    0.001   33.344    1.042 /lma1/mpound/dysh/src/dysh/util/selection.py:1281(read)
      448    0.017    0.000   27.378    0.061 /lma1/mpound/dysh/src/dysh/util/selection.py:1124(flag)
      448    0.009    0.000   23.104    0.052 /lma1/mpound/dysh/src/dysh/util/selection.py:534(_base_select)
      448    0.014    0.000   12.623    0.028 /lma1/mpound/dysh/src/dysh/util/selection.py:458(_addrow)
      928    0.140    0.000   11.608    0.013 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3575(sort)
   403835    1.225    0.000   11.171    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:704(__array_finalize__)
      448    0.022    0.000   10.392    0.023 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
      448    0.005    0.000    9.655    0.022 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
      448    0.339    0.001    8.524    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
    88608    0.324    0.000    7.842    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1362(__setitem__)
   216272    0.307    0.000    7.130    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:722(__array_wrap__)
4336961/3070254    4.333    0.000    6.994    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1270(__setattr__)
      452    0.001    0.000    6.887    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3153(add_row)
 1348/452    0.137    0.000    6.885    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3206(insert_row)
    86752    0.439    0.000    6.879    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1338(_check_string_truncate)
   216272    0.190    0.000    6.693    0.000 {function BaseColumn.__array_wrap__ at 0x7f544b37fba0}
   403835    0.692    0.000    6.299    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1124(_copy_attrs)
    94304    0.793    0.000    6.000    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
    44408    0.161    0.000    5.569    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1384(insert)
       32    0.001    0.000    5.123    0.160 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
2014457/1863492    0.666    0.000    5.105    0.000 {built-in method builtins.setattr}
    43904    0.109    0.000    4.287    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
       32    0.003    0.000    4.072    0.127 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
      448    0.001    0.000    4.051    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:341(_check_numbers)
42112/448    0.476    0.000    4.050    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:378(_check_type)
   184570    0.682    0.000    3.306    0.000 {method 'reduce' of 'numpy.ufunc' objects}
501374/494778    0.257    0.000    3.263    0.000 {method 'view' of 'numpy.ndarray' objects}
   129464    0.051    0.000    3.195    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/_methods.py:42(_amax)
    86900    0.055    0.000    3.093    0.000 {method 'max' of 'numpy.ndarray' objects}
    45168    0.514    0.000    2.870    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/_function_base_impl.py:5456(insert)
   133280    0.286    0.000    2.812    0.000 {method 'take' of 'numpy.ndarray' objects}
6115006/6114926    1.261    0.000    2.333    0.000 {built-in method builtins.getattr}
      136    0.004    0.000    2.296    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
10381346/10371542    1.779    0.000    2.117    0.000 {built-in method builtins.isinstance}
    45168    0.077    0.000    1.710    0.000 {method 'wrap' of 'numpy._core._multiarray_umath._array_converter' objects}
47396/23364    0.054    0.000    1.616    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/decorators.py:842(__get__)
    10064    0.166    0.000    1.606    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
       32    0.000    0.000    1.549    0.048 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:404(data)
       32    0.000    0.000    1.528    0.048 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:174(_get_tbdata)
     2408    0.008    0.000    1.509    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
    98224    0.437    0.000    1.390    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
    42148    0.229    0.000    1.287    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/ma/core.py:2879(__new__)
       68    0.004    0.000    1.261    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
       32    0.000    0.000    1.141    0.036 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       32    0.002    0.000    1.140    0.036 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       32    0.000    0.000    1.139    0.036 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
    83328    0.241    0.000    1.137    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/ma/core.py:3293(__getitem__)
