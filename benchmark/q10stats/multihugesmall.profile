Loading not more than 5 from /lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT14B_480_06/AGBT14B_480_06.raw.vegas
Will load 5 of 5 files. FITS size per file 3713.27MB, Flag lines 6
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 5 FITS files
Found ostype= linux
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 5 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 5 FITS files
Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.
Loaded 5 FITS files
# Dysh Benchmark: GBTFITSLoad
 name  time   VmSize VmRSS  #files file_size totsize  nchan   nrow nIF nFd nPol #flags skipflags
        ms    Mbyte  Mbyte           Mbyte    Mbyte                                             
----- ------ ------- ------ ------ --------- -------- ------ ----- --- --- ---- ------ ---------
load1 4493.4 22916.0 7032.2      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load2 4148.5 23094.4 7210.2      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load3 4005.0 23074.1 7190.5      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load4 4057.8 23074.3 7190.7      5   3713.27 18566.35 131072 35370   5   1    2      6     False
         20691090 function calls (20488971 primitive calls) in 16.617 seconds

   Ordered by: cumulative time
   List reduced from 1646 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   16.124    4.031 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.002    0.000   16.122    4.030 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
       20    0.001    0.000   10.613    0.531 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
       20    0.013    0.001    8.182    0.409 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
        4    0.059    0.015    5.429    1.357 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:892(_create_index_if_needed)
        4    0.152    0.038    4.061    1.015 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:970(_construct_integration_number)
      944    0.990    0.001    3.977    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       20    0.000    0.000    3.751    0.188 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       20    0.007    0.000    3.750    0.187 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       20    0.000    0.000    3.743    0.187 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       20    0.002    0.000    3.706    0.185 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      460    0.002    0.000    3.668    0.008 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
    20156    0.071    0.000    3.606    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
    20156    0.022    0.000    3.525    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1149(__iter__)
    19952    1.096    0.000    3.069    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1180(_chop)
       20    0.000    0.000    2.890    0.145 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:505(drop_fields)
       20    0.878    0.044    2.888    0.144 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:32(recursive_fill_fields)
     1508    0.006    0.000    2.653    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      484    0.007    0.000    2.556    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
       20    0.001    0.000    2.333    0.117 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:75(_init_flags)
     1484    0.008    0.000    2.007    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:692(field)
     1484    0.007    0.000    1.950    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:883(_convert_other)
      460    0.002    0.000    1.924    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/util.py:254(decode_ascii)
      460    0.021    0.000    1.921    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:509(decode)
      460    0.004    0.000    1.843    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
     2424    1.794    0.001    1.800    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/numeric.py:290(full)
      460    0.001    0.000    1.724    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
      460    0.001    0.000    1.673    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
       88    0.003    0.000    1.433    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
    19952    1.375    0.000    1.420    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexes/base.py:5425(_getitem_slice)
      460    1.392    0.003    1.392    0.003 {built-in method numpy._core._multiarray_umath._vec_string}
  3254040    0.640    0.000    1.189    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
  3254040    0.634    0.000    1.081    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
     6512    0.107    0.000    1.048    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
30296/15276    0.034    0.000    0.965    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/decorators.py:842(__get__)
       20    0.000    0.000    0.920    0.046 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:404(data)
       20    0.000    0.000    0.906    0.045 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:174(_get_tbdata)
        4    0.003    0.001    0.873    0.218 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:944(_construct_procedure)
       24    0.000    0.000    0.824    0.034 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:834(split)
13710/9442    0.010    0.000    0.803    0.000 {method 'view' of 'numpy.ndarray' objects}
       44    0.003    0.000    0.775    0.018 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
       44    0.020    0.000    0.722    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
       24    0.000    0.000    0.648    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:555(nchan)
       24    0.000    0.000    0.648    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:415(rawspectrum)
214361/116676    0.060    0.000    0.603    0.000 {built-in method builtins.setattr}
  3254056    0.549    0.000    0.549    0.000 {method 'decode' of 'bytes' objects}
   104192    0.104    0.000    0.545    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
       24    0.000    0.000    0.514    0.021 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:327(_str_split)
       20    0.000    0.000    0.508    0.025 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:395(columns)
      460    0.007    0.000    0.507    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:98(_to_bytes_or_str_array)


