Using file f1=PosixPath('/lma1/teuben/GBT/dysh_data/acceptance_testing/data/AGBT20B_336_01/AGBT20B_336_01.raw.vegas')
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
 name   time  VmSize VmRSS  #files file_size totsize nchan  nrow nIF nFd nPol #flags skipflags
         ms   Mbyte  Mbyte           Mbyte    Mbyte                                           
----- ------- ------ ------ ------ --------- ------- ----- ----- --- --- ---- ------ ---------
load1 15768.1 4224.2 1335.0      8     56.18  449.44  1024 95136   4   2    2      7     False
load2 15822.2 4673.5 1786.3      8     56.18  449.44  1024 95136   4   2    2      7     False
load3 15748.9 4755.4 1867.3      8     56.18  449.44  1024 95136   4   2    2      7     False
load4 15696.6 4799.2 1910.6      8     56.18  449.44  1024 95136   4   2    2      7     False
         116004696 function calls (115544104 primitive calls) in 62.688 seconds

   Ordered by: cumulative time
   List reduced from 1723 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   61.913   15.478 /lma1/mpound/dysh/src/dysh/log.py:337(wrapper)
        4    0.001    0.000   61.910   15.478 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:69(__init__)
        4    0.237    0.059   41.264   10.316 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:861(_create_index_if_needed)
        8    0.000    0.000   33.462    4.183 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
        8    0.000    0.000   33.360    4.170 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
 1520/784    0.009    0.000   33.049    0.042 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
     1488    0.004    0.000   32.990    0.022 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
    12020    0.180    0.000   32.788    0.003 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
     1488    0.006    0.000   32.702    0.022 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   761096    0.251    0.000   28.312    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:263(self_iter)
   761088    0.667    0.000   28.061    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/shapes.py:243(__getitem__)
   761088    8.359    0.000   27.394    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/core.py:1321(_apply)
       32    0.001    0.000   20.469    0.640 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
       32    0.033    0.001   19.387    0.606 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
  1522208    2.488    0.000   11.796    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:156(__init__)
     1508    2.640    0.002   10.272    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       32    0.000    0.000    9.767    0.305 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       32    0.025    0.001    9.766    0.305 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       32    0.000    0.000    9.740    0.304 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       32    0.003    0.000    9.662    0.302 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      736    0.003    0.000    9.597    0.013 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
       32    0.000    0.000    6.850    0.214 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:505(drop_fields)
       32    2.240    0.070    6.847    0.214 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/lib/recfunctions.py:32(recursive_fill_fields)
      772    0.019    0.000    6.429    0.008 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
     2408    0.008    0.000    5.592    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      736    0.008    0.000    4.850    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
     2372    0.013    0.000    4.600    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:692(field)
      736    0.002    0.000    4.511    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
     2372    0.010    0.000    4.511    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:883(_convert_other)
      736    0.003    0.000    4.468    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/util.py:254(decode_ascii)
      736    0.053    0.000    4.464    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/numpy/_core/strings.py:509(decode)
      736    0.001    0.000    4.426    0.006 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
       72    0.002    0.000    4.157    0.058 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:123(maybe_convert_platform)
       72    1.098    0.015    4.155    0.058 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:1579(construct_1d_object_array_from_listlike)
  3805488    3.899    0.000    3.899    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/core.py:1758(__getattr__)
  3044368    1.912    0.000    3.636    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:2503(_validate_jd_for_storage)
  1522184    1.136    0.000    3.182    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:242(jd2)
  8752512    1.734    0.000    3.171    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
3915911/3915839    1.171    0.000    3.165    0.000 {built-in method builtins.getattr}
      736    3.082    0.004    3.082    0.004 {built-in method numpy._core._multiarray_umath._vec_string}
       32    0.003    0.000    2.949    0.092 /lma1/mpound/dysh/src/dysh/util/selection.py:1211(read)
  8752512    1.684    0.000    2.852    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
  1522184    0.620    0.000    2.544    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/time/formats.py:232(jd1)
      136    0.004    0.000    2.276    0.017 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
       32    0.001    0.000    2.194    0.069 /lma1/mpound/dysh/src/dysh/util/selection.py:1060(flag)
        4    0.061    0.015    2.104    0.526 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:939(_construct_integration_number)
        4    0.010    0.003    2.009    0.502 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:913(_construct_procedure)
14618884/14611872    1.940    0.000    2.000    0.000 {built-in method builtins.isinstance}
       36    0.000    0.000    1.898    0.053 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/strings/accessor.py:834(split)
    11904    0.038    0.000    1.894    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
