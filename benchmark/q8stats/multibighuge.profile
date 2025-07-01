Using file f1=PosixPath('/lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT17B_319_06/AGBT17B_319_06.raw.vegas')
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
load1 37816.6 5962.9 1980.2      3    768.09 2304.27 32768 17496   3   2    2    102     False
load2 35816.9 6077.0 2093.5      3    768.09 2304.27 32768 17496   3   2    2    102     False
load3 35391.4 6078.4 2094.8      3    768.09 2304.27 32768 17496   3   2    2    102     False
load4 35572.2 6091.0 2107.2      3    768.09 2304.27 32768 17496   3   2    2    102     False
         264116031 function calls (258046255 primitive calls) in 144.458 seconds

   Ordered by: cumulative time
   List reduced from 1739 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000  144.090   36.022 /lma1/mpound/dysh/src/dysh/log.py:337(wrapper)
        4    0.002    0.000  144.086   36.022 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:69(__init__)
        4    0.021    0.005  136.779   34.195 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:861(_create_index_if_needed)
       12    0.038    0.003  128.588   10.716 /lma1/mpound/dysh/src/dysh/util/selection.py:1211(read)
     1152    0.023    0.000  113.178    0.098 /lma1/mpound/dysh/src/dysh/util/selection.py:1060(flag)
     1152    0.021    0.000  111.125    0.096 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
     1152    0.012    0.000   84.201    0.073 /lma1/mpound/dysh/src/dysh/util/selection.py:446(_addrow)
     1152    0.350    0.000   50.779    0.044 /lma1/mpound/dysh/src/dysh/util/selection.py:419(_check_for_duplicates)
  1348551    4.104    0.000   37.120    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:704(__array_finalize__)
     2316    0.350    0.000   30.914    0.013 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3575(sort)
   174572    8.265    0.000   29.598    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:1024(indices)
   166464    0.119    0.000   28.793    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:1037(loc)
   166464    0.093    0.000   28.675    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:825(__init__)
     1152    0.055    0.000   26.704    0.023 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
     1152    0.013    0.000   24.819    0.022 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
14060273/9877250   14.233    0.000   23.162    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1270(__setattr__)
     1152    0.858    0.001   22.119    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
  1348551    2.359    0.000   21.324    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1124(_copy_attrs)
   221172    0.897    0.000   21.261    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1362(__setitem__)
     1156    0.011    0.000   17.809    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3153(add_row)
3460/1156    0.449    0.000   17.797    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3206(insert_row)
   544136    0.778    0.000   17.797    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:722(__array_wrap__)
   216540    1.373    0.000   17.633    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1338(_check_string_truncate)
   166464    0.254    0.000   16.762    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:863(__getitem__)
   544136    0.482    0.000   16.681    0.000 {function BaseColumn.__array_wrap__ at 0x7f9387bc2160}
   166464    0.585    0.000   15.715    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:829(_get_rows)
   239196    2.010    0.000   15.035    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
5886249/5824084    2.019    0.000   14.915    0.000 {built-in method builtins.setattr}
   114104    0.396    0.000   14.109    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1384(insert)
   166464    0.129    0.000   13.171    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:515(find)
   166464    0.168    0.000   13.013    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:250(find)
 18663463    6.502    0.000   12.982    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:216(__get__)
   166464    0.848    0.000   12.844    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/sorted_array.py:130(find)
 18602569    8.722    0.000   11.686    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/data_info.py:340(__get__)
   112928    0.733    0.000   11.251    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
   240844    1.109    0.000    9.518    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
32332295/32332223    5.795    0.000    9.322    0.000 {built-in method builtins.getattr}
   519390    2.004    0.000    8.626    0.000 {method 'reduce' of 'numpy.ufunc' objects}
   325616    0.121    0.000    8.120    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/_methods.py:42(_amax)
   216172    0.136    0.000    7.840    0.000 {method 'max' of 'numpy.ndarray' objects}
   339376    0.879    0.000    7.254    0.000 {method 'take' of 'numpy.ndarray' objects}
       12    0.000    0.000    7.209    0.601 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
   114424    1.364    0.000    7.190    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/_function_base_impl.py:5456(insert)
       12    0.006    0.000    6.599    0.550 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
        8    0.000    0.000    6.420    0.803 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
  600/324    0.003    0.000    6.201    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
        8    0.000    0.000    6.196    0.775 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
      588    0.002    0.000    6.179    0.011 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
31403813/31394629    5.327    0.000    6.150    0.000 {built-in method builtins.isinstance}
      588    0.002    0.000    6.068    0.010 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
