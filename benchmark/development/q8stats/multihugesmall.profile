Using file f1=PosixPath('/lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT14B_480_06/AGBT14B_480_06.raw.vegas')
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
 name   time   VmSize VmRSS  #files file_size totsize  nchan   nrow nIF nFd nPol #flags skipflags
         ms    Mbyte  Mbyte           Mbyte    Mbyte                                             
----- ------- ------- ------ ------ --------- -------- ------ ----- --- --- ---- ------ ---------
load1 20655.3 25428.8 7068.8      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load2  6939.9 25597.3 7237.4      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load3  7017.5 25600.6 7240.6      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load4  7057.7 25601.0 7241.0      5   3713.27 18566.35 131072 35370   5   1    2      6     False
         44491906 function calls (44290784 primitive calls) in 41.552 seconds

   Ordered by: cumulative time
   List reduced from 1542 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   41.062   10.265 /lma1/mpound/dysh/src/dysh/log.py:337(wrapper)
        4    0.002    0.000   41.060   10.265 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:69(__init__)
       20    0.001    0.000   24.208    1.210 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
       20    0.011    0.001   21.814    1.091 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
        4    0.052    0.013   16.774    4.193 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:861(_create_index_if_needed)
       20    0.000    0.000   16.550    0.828 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:505(drop_fields)
       20    0.881    0.044   16.548    0.827 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:32(recursive_fill_fields)
     1508    0.005    0.000   16.319    0.011 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
     1484    0.008    0.000   15.663    0.011 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:692(field)
     1484    0.007    0.000   15.607    0.011 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:883(_convert_other)
      460    0.002    0.000   15.580    0.034 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/util.py:254(decode_ascii)
      460    0.020    0.000   15.577    0.034 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:509(decode)
      460   15.060    0.033   15.060    0.033 {built-in method numpy._core._multiarray_umath._vec_string}
        8    0.000    0.000   12.595    1.574 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
        8    0.000    0.000   12.496    1.562 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
  968/508    0.006    0.000   12.455    0.025 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
      948    0.003    0.000   12.420    0.013 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
     3612    0.063    0.000   12.270    0.003 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
      948    0.003    0.000   12.247    0.013 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   282968    0.092    0.000   10.614    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:263(self_iter)
   282960    0.236    0.000   10.522    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:243(__getitem__)
   282960    3.225    0.000   10.285    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/core.py:1321(_apply)
   565952    0.942    0.000    4.382    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:156(__init__)
       20    0.000    0.000    3.742    0.187 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       20    0.006    0.000    3.742    0.187 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
      944    1.000    0.001    3.736    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       20    0.000    0.000    3.735    0.187 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       20    0.002    0.000    3.698    0.185 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      460    0.002    0.000    3.661    0.008 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
        4    0.141    0.035    3.244    0.811 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:939(_construct_integration_number)
    20156    0.064    0.000    2.808    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
    20156    0.020    0.000    2.733    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1149(__iter__)
      484    0.006    0.000    2.302    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
       20    0.001    0.000    2.275    0.114 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:75(_init_flags)
    19952    1.064    0.000    2.201    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1180(_chop)
      460    0.003    0.000    1.845    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
     2424    1.729    0.001    1.734    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/numeric.py:290(full)
      460    0.001    0.000    1.715    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
      460    0.001    0.000    1.666    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
       48    0.001    0.000    1.548    0.032 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:123(maybe_convert_platform)
       48    0.401    0.008    1.546    0.032 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:1579(construct_1d_object_array_from_listlike)
  1414848    1.459    0.000    1.459    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/core.py:1758(__getattr__)
       88    0.003    0.000    1.458    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
  1131856    0.717    0.000    1.352    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:2503(_validate_jd_for_storage)
  3254040    0.639    0.000    1.185    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
   565928    0.416    0.000    1.170    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:242(jd2)
1420667/1420595    0.419    0.000    1.120    0.000 {built-in method builtins.getattr}
  3254040    0.623    0.000    1.069    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
     6512    0.112    0.000    1.067    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
    19952    0.949    0.000    0.992    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/indexes/base.py:5425(_getitem_slice)
