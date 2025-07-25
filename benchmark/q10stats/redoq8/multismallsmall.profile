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
load1 15731.7 2015.2 1398.0      8     56.18  449.44  1024 95136   4   2    2      7     False
load2 15486.8 2456.7 1840.3      8     56.18  449.44  1024 95136   4   2    2      7     False
load3 15621.6 2543.6 1926.6      8     56.18  449.44  1024 95136   4   2    2      7     False
load4 15487.2 2581.3 1963.4      8     56.18  449.44  1024 95136   4   2    2      7     False
         115002746 function calls (114617908 primitive calls) in 61.949 seconds

   Ordered by: cumulative time
   List reduced from 1885 to 100 due to restriction <100>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   61.215   15.304 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.001    0.000   61.212   15.303 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
        4    0.235    0.059   40.821   10.205 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:870(_create_index_if_needed)
        8    0.000    0.000   33.734    4.217 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
        8    0.000    0.000   33.633    4.204 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
 1520/784    0.009    0.000   33.306    0.042 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
     1488    0.004    0.000   33.245    0.022 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
    11892    0.168    0.000   33.046    0.003 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
     1488    0.006    0.000   32.958    0.022 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   761096    0.250    0.000   28.604    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:252(self_iter)
   761088    0.630    0.000   28.353    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:232(__getitem__)
   761088    8.671    0.000   27.723    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
       32    0.001    0.000   20.235    0.632 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
       32    0.034    0.001   19.166    0.599 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
  1522208    2.559    0.000   11.717    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
     1508    2.588    0.002    9.857    0.007 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       32    0.001    0.000    9.591    0.300 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       32    0.026    0.001    9.589    0.300 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       32    0.000    0.000    9.563    0.299 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       32    0.003    0.000    9.484    0.296 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      736    0.003    0.000    9.417    0.013 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
       32    0.001    0.000    6.673    0.209 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:501(drop_fields)
       32    2.057    0.064    6.670    0.208 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:35(recursive_fill_fields)
      772    0.019    0.000    6.084    0.008 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
     2408    0.008    0.000    5.575    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:503(__getitem__)
      736    0.012    0.000    4.787    0.007 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
     2372    0.013    0.000    4.607    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:684(field)
     2372    0.011    0.000    4.518    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:875(_convert_other)
      736    0.003    0.000    4.475    0.006 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/util.py:253(decode_ascii)
      736    0.051    0.000    4.470    0.006 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:572(decode)
      736    0.002    0.000    4.427    0.006 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
      736    0.001    0.000    4.347    0.006 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
      136    0.003    0.000    4.139    0.030 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:124(maybe_convert_platform)
      136    1.119    0.008    4.135    0.030 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:1580(construct_1d_object_array_from_listlike)
  3805488    3.887    0.000    3.887    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1792(__getattr__)
  3044368    1.956    0.000    3.631    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:2453(_validate_jd_for_storage)
3778962/3778890    1.260    0.000    3.171    0.000 {built-in method builtins.getattr}
  8752512    1.696    0.000    3.170    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
  1522184    1.103    0.000    3.114    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:241(jd2)
    12336    3.025    0.000    3.025    0.000 {built-in method numpy.core._multiarray_umath._vec_string}
  8752512    1.660    0.000    2.842    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
       32    0.003    0.000    2.607    0.081 /lma1/mpound/dysh/src/dysh/util/selection.py:1213(read)
  1522184    0.654    0.000    2.601    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:231(jd1)
      136    0.004    0.000    2.373    0.017 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1466(__init__)
       32    0.001    0.000    2.238    0.070 /lma1/mpound/dysh/src/dysh/util/selection.py:1062(flag)
        4    0.068    0.017    1.984    0.496 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:948(_construct_integration_number)
14399145/14391973    1.876    0.000    1.938    0.000 {built-in method builtins.isinstance}
       32    0.001    0.000    1.929    0.060 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
    11904    0.036    0.000    1.771    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
        4    0.010    0.003    1.746    0.436 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:922(_construct_procedure)
    11904    0.011    0.000    1.724    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1149(__iter__)
47396/23364    0.054    0.000    1.683    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/decorators.py:827(__get__)
  4566600    1.198    0.000    1.669    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:441(_select_subfmts)
       36    0.000    0.000    1.637    0.045 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:834(split)
       32    0.000    0.000    1.614    0.050 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:405(data)
    10064    0.162    0.000    1.607    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
       32    0.000    0.000    1.593    0.050 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:175(_get_tbdata)
47330/40734    0.028    0.000    1.505    0.000 {method 'view' of 'numpy.ndarray' objects}
  8752527    1.474    0.000    1.474    0.000 {method 'decode' of 'bytes' objects}
      736    0.014    0.000    1.422    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:60(_to_bytes_or_str_array)
32192/30616    1.390    0.000    1.393    0.000 {built-in method numpy.asarray}
       68    0.004    0.000    1.335    0.020 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:230(__array_finalize__)
    11896    0.647    0.000    1.250    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1180(_chop)
1162521/1011556    0.253    0.000    1.236    0.000 {built-in method builtins.setattr}
  1545906    1.231    0.000    1.231    0.000 {built-in method numpy.array}
  1522208    0.642    0.000    1.210    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:211(in_subfmt)
  8815192    1.196    0.000    1.196    0.000 {method 'strip' of 'str' objects}
       32    0.002    0.000    1.196    0.037 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
  1522176    0.567    0.000    1.145    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:193(_get_allowed_subfmt)
       68    0.029    0.000    1.102    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1524(_init_from_array)
  1522208    0.569    0.000    1.091    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:221(out_subfmt)
   761112    0.492    0.000    0.997    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:223(__len__)
       32    0.000    0.000    0.976    0.031 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
       36    0.000    0.000    0.971    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:544(nchan)
       36    0.000    0.000    0.970    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:404(rawspectrum)
     1508    0.304    0.000    0.966    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:255(_wrap_result)
       32    0.001    0.000    0.919    0.029 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:75(_init_flags)
       32    0.025    0.001    0.851    0.027 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
   161024    0.152    0.000    0.844    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
       32    0.000    0.000    0.822    0.026 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:396(columns)
       36    0.000    0.000    0.811    0.023 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:327(_str_split)
       32    0.025    0.001    0.808    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1580(_init_from_table)
  1522208    0.611    0.000    0.785    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:276(precision)
        8    0.257    0.032    0.769    0.096 /lma1/mpound/dysh/src/dysh/util/core.py:103(gbt_timestamp_to_time)
      148    0.003    0.000    0.761    0.005 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:694(__init__)
   761096    0.246    0.000    0.741    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:777(scale)
    63310    0.039    0.000    0.725    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:333(__get__)
    39344    0.049    0.000    0.724    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/card.py:282(value)
       32    0.000    0.000    0.723    0.023 /lma1/mpound/dysh/src/dysh/util/selection.py:446(_addrow)
        4    0.000    0.000    0.691    0.173 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:314(stats)
    63310    0.159    0.000    0.686    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:442(__call__)
   761088    0.472    0.000    0.676    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:358(<lambda>)
       40    0.001    0.000    0.675    0.017 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
      100    0.022    0.000    0.669    0.007 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/construction.py:96(arrays_to_mgr)
       40    0.000    0.000    0.663    0.017 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:869(take)
       40    0.003    0.000    0.656    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:623(reindex_indexer)
     2792    0.005    0.000    0.650    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/blocks.py:1287(take_nd)
     2800    0.004    0.000    0.642    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:59(take_nd)
     2800    0.522    0.000    0.637    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:120(_take_nd_ndarray)
    10064    0.020    0.000    0.607    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:871(name)


