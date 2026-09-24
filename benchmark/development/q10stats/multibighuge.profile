Loading not more than 3 from /lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT17B_319_06/AGBT17B_319_06.raw.vegas
Will load 3 of 3 files. FITS size per file 768.09MB, Flag lines 102
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 3 FITS files
Found ostype= linux
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 3 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 3 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 3 FITS files
# Dysh Benchmark: GBTFITSLoad
 name   time  VmSize VmRSS  #files file_size totsize nchan  nrow nIF nFd nPol #flags skipflags
         ms   Mbyte  Mbyte           Mbyte    Mbyte                                           
----- ------- ------ ------ ------ --------- ------- ----- ----- --- --- ---- ------ ---------
load1 21384.4 3467.5 1963.0      3    768.09 2304.27 32768 17496   3   2    2    102     False
load2 21325.0 3577.0 2072.2      3    768.09 2304.27 32768 17496   3   2    2    102     False
load3 21448.4 3567.9 2064.2      3    768.09 2304.27 32768 17496   3   2    2    102     False
load4 21356.3 3569.6 2065.8      3    768.09 2304.27 32768 17496   3   2    2    102     False
         136112519 function calls (131527258 primitive calls) in 85.443 seconds

   Ordered by: cumulative time
   List reduced from 1835 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   85.139   21.285 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.001    0.000   85.137   21.284 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
        4    0.015    0.004   80.180   20.045 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:892(_create_index_if_needed)
       12    0.039    0.003   78.389    6.532 /lma1/mpound/dysh/src/dysh/util/selection.py:1281(read)
     1152    0.023    0.000   62.757    0.054 /lma1/mpound/dysh/src/dysh/util/selection.py:1124(flag)
     1152    0.024    0.000   60.716    0.053 /lma1/mpound/dysh/src/dysh/util/selection.py:534(_base_select)
     1152    0.036    0.000   33.864    0.029 /lma1/mpound/dysh/src/dysh/util/selection.py:458(_addrow)
     2316    0.365    0.000   31.374    0.014 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3575(sort)
  1017927    3.127    0.000   28.647    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:704(__array_finalize__)
     1152    0.056    0.000   26.622    0.023 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
     1152    0.012    0.000   24.714    0.021 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
     1152    0.858    0.001   22.030    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
   221172    0.881    0.000   21.578    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1362(__setitem__)
   544136    0.796    0.000   18.172    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:722(__array_wrap__)
   216540    1.412    0.000   17.957    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1338(_check_string_truncate)
10919345/7728194   11.231    0.000   17.957    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1270(__setattr__)
     1156    0.003    0.000   17.907    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3153(add_row)
3460/1156    0.355    0.000   17.902    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3206(insert_row)
   544136    0.489    0.000   17.041    0.000 {function BaseColumn.__array_wrap__ at 0x7fa45a4dfba0}
  1017927    1.753    0.000   16.349    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1124(_copy_attrs)
   239228    2.003    0.000   14.979    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
   114104    0.416    0.000   14.440    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1384(insert)
4423881/4361716    1.503    0.000   11.347    0.000 {built-in method builtins.setattr}
   112928    0.710    0.000   11.166    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
   369998    1.480    0.000    8.218    0.000 {method 'reduce' of 'numpy.ufunc' objects}
   326972    0.121    0.000    8.211    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/_methods.py:42(_amax)
   217528    0.136    0.000    7.939    0.000 {method 'max' of 'numpy.ndarray' objects}
   339376    0.882    0.000    7.377    0.000 {method 'take' of 'numpy.ndarray' objects}
   114424    1.361    0.000    7.299    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/_function_base_impl.py:5456(insert)
12781326/12781246    2.877    0.000    5.692    0.000 {built-in method builtins.getattr}
       12    0.000    0.000    4.864    0.405 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
23131223/23120831    3.923    0.000    4.590    0.000 {built-in method builtins.isinstance}
       12    0.005    0.000    4.265    0.355 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
   114424    0.197    0.000    4.227    0.000 {method 'wrap' of 'numpy._core._multiarray_umath._array_converter' objects}
776990/774274    0.390    0.000    3.593    0.000 {method 'view' of 'numpy.ndarray' objects}
   240884    1.075    0.000    3.428    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
   886405    0.659    0.000    2.644    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:572(__set__)
   112928    0.919    0.000    2.482    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:287(get_dtypes)
   903051    0.490    0.000    2.313    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/metadata/core.py:91(__get__)
   221172    0.896    0.000    2.262    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:607(adjust_indices)
  2793511    1.004    0.000    2.258    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:216(__get__)
  1017927    0.799    0.000    2.204    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:764(name)
     1152    0.013    0.000    2.170    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/computation/eval.py:170(eval)
     1348    0.009    0.000    2.050    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
     2304    0.004    0.000    2.043    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:540(insert_row)
     2304    0.210    0.000    2.038    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:164(insert_row)
   239160    0.552    0.000    2.016    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:1863(from_array)
  2567305    1.519    0.000    2.004    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:340(__get__)
   671389    1.246    0.000    1.985    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:355(__set__)
     1184    0.004    0.000    1.851    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexing.py:1176(__getitem__)
