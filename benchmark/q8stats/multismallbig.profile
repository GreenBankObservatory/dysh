Using file f1=PosixPath('/lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT23A_432_03/AGBT23A_432_03.raw.vegas')
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
 name   time  VmSize VmRSS  #files file_size totsize nchan nrow nIF nFd nPol #flags skipflags
         ms   Mbyte  Mbyte           Mbyte    Mbyte                                          
----- ------- ------ ------ ------ --------- ------- ----- ---- --- --- ---- ------ ---------
load1 12756.8 3790.8  904.3      8     66.72  533.76 16384 8064   4   2    2     20     False
load2 12729.1 3855.8  968.7      8     66.72  533.76 16384 8064   4   2    2     20     False
load3 12585.1 3854.9  967.7      8     66.72  533.76 16384 8064   4   2    2     20     False
load4 12553.0 3973.8 1086.8      8     66.72  533.76 16384 8064   4   2    2     20     False
         87102780 function calls (84655628 primitive calls) in 50.583 seconds

   Ordered by: cumulative time
   List reduced from 1723 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   50.364   12.591 /lma1/mpound/dysh/src/dysh/log.py:337(wrapper)
        4    0.002    0.000   50.361   12.590 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:69(__init__)
        4    0.015    0.004   45.148   11.287 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:861(_create_index_if_needed)
       32    0.037    0.001   41.448    1.295 /lma1/mpound/dysh/src/dysh/util/selection.py:1211(read)
      448    0.017    0.000   35.352    0.079 /lma1/mpound/dysh/src/dysh/util/selection.py:1060(flag)
      448    0.008    0.000   31.056    0.069 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
      448    0.004    0.000   20.654    0.046 /lma1/mpound/dysh/src/dysh/util/selection.py:446(_addrow)
   453563    1.388    0.000   12.859    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:704(__array_finalize__)
      928    0.146    0.000   11.855    0.013 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3575(sort)
      448    0.021    0.000   10.316    0.023 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
      448    0.005    0.000    9.599    0.021 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
      448    0.335    0.001    8.444    0.019 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
4809377/3393486    5.049    0.000    8.087    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1270(__setattr__)
    88608    0.323    0.000    8.007    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1362(__setitem__)
      448    0.054    0.000    7.769    0.017 /lma1/mpound/dysh/src/dysh/util/selection.py:419(_check_for_duplicates)
   216272    0.316    0.000    7.344    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:722(__array_wrap__)
   453563    0.788    0.000    7.297    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1124(_copy_attrs)
      452    0.004    0.000    7.072    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3153(add_row)
 1348/452    0.175    0.000    7.067    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:3206(insert_row)
    86752    0.447    0.000    7.044    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1338(_check_string_truncate)
   216272    0.193    0.000    6.893    0.000 {function BaseColumn.__array_wrap__ at 0x7f775a1c2160}
    94272    0.790    0.000    5.953    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
2277785/2126820    0.768    0.000    5.814    0.000 {built-in method builtins.setattr}
    44408    0.154    0.000    5.639    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/column.py:1384(insert)
       32    0.001    0.000    5.176    0.162 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
    28532    1.316    0.000    4.778    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:1024(indices)
    25312    0.019    0.000    4.402    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/table.py:1037(loc)
    25312    0.015    0.000    4.383    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:825(__init__)
    43904    0.118    0.000    4.190    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
    98184    0.448    0.000    4.188    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
       32    0.003    0.000    4.119    0.129 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
      448    0.001    0.000    4.072    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:333(_check_numbers)
42112/448    0.485    0.000    4.071    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:366(_check_type)
   208442    0.774    0.000    3.468    0.000 {method 'reduce' of 'numpy.ufunc' objects}
   128948    0.048    0.000    3.274    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/_methods.py:42(_amax)
550086/543490    0.282    0.000    3.241    0.000 {method 'view' of 'numpy.ndarray' objects}
    86384    0.056    0.000    3.170    0.000 {method 'max' of 'numpy.ndarray' objects}
9230839/9230767    1.781    0.000    3.132    0.000 {built-in method builtins.getattr}
 1520/784    0.008    0.000    3.097    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
     1488    0.004    0.000    3.046    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
        8    0.000    0.000    2.955    0.369 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
    45168    0.515    0.000    2.934    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/_function_base_impl.py:5456(insert)
   133280    0.299    0.000    2.880    0.000 {method 'take' of 'numpy.ndarray' objects}
        8    0.000    0.000    2.855    0.357 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
     1488    0.004    0.000    2.815    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
    25312    0.039    0.000    2.613    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:863(__getitem__)
12416160/12406860    2.126    0.000    2.502    0.000 {built-in method builtins.isinstance}
    25312    0.096    0.000    2.452    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/index.py:829(_get_rows)
    64520    0.022    0.000    2.414    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:263(self_iter)
    64512    0.054    0.000    2.391    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:243(__getitem__)
