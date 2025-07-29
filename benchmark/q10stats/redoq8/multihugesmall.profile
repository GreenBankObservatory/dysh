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
load1  7787.2 23222.0 7131.2      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load2 10998.6 23387.9 7296.8      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load3  6928.3 23383.2 7292.3      5   3713.27 18566.35 131072 35370   5   1    2      6     False
load4  9502.7 23400.2 7309.0      5   3713.27 18566.35 131072 35370   5   1    2      6     False
         44491366 function calls (44314398 primitive calls) in 35.094 seconds

   Ordered by: cumulative time
   List reduced from 1612 to 100 due to restriction <100>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   34.593    8.648 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.001    0.000   34.591    8.648 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
       20    0.001    0.000   17.737    0.887 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
        4    0.053    0.013   16.783    4.196 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:870(_create_index_if_needed)
        8    0.000    0.000   12.636    1.580 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
        8    0.000    0.000   12.538    1.567 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
  968/508    0.006    0.000   12.486    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
      948    0.003    0.000   12.451    0.013 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
     3612    0.061    0.000   12.306    0.003 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
      948    0.003    0.000   12.282    0.013 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   282968    0.097    0.000   10.624    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:252(self_iter)
   282960    0.235    0.000   10.527    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:232(__getitem__)
   282960    3.281    0.000   10.293    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
       20    0.001    0.000    9.744    0.487 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:75(_init_flags)
     2424    9.211    0.004    9.217    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/numeric.py:274(full)
       20    0.012    0.001    7.894    0.395 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
   565952    0.917    0.000    4.304    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
      944    0.958    0.001    3.684    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       20    0.000    0.000    3.653    0.183 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       20    0.007    0.000    3.653    0.183 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       20    0.000    0.000    3.646    0.182 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       20    0.002    0.000    3.610    0.181 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      460    0.002    0.000    3.573    0.008 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
        4    0.152    0.038    3.196    0.799 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:948(_construct_integration_number)
    20156    0.062    0.000    2.754    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:607(get_iterator)
       20    0.000    0.000    2.747    0.137 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:501(drop_fields)
       20    0.799    0.040    2.745    0.137 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:35(recursive_fill_fields)
    20156    0.020    0.000    2.681    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1149(__iter__)
     1508    0.005    0.000    2.582    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:503(__getitem__)
      484    0.007    0.000    2.280    0.005 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:129(wrapper)
    19952    1.028    0.000    2.167    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1180(_chop)
     1484    0.008    0.000    1.942    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:684(field)
     1484    0.007    0.000    1.887    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:875(_convert_other)
      460    0.002    0.000    1.860    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/util.py:253(decode_ascii)
      460    0.020    0.000    1.857    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:572(decode)
      460    0.003    0.000    1.810    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:1972(decode)
      460    0.001    0.000    1.677    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2115(strip)
      460    0.001    0.000    1.629    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:450(_str_strip)
       48    0.001    0.000    1.575    0.033 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:124(maybe_convert_platform)
       48    0.419    0.009    1.574    0.033 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:1580(construct_1d_object_array_from_listlike)
  1414848    1.448    0.000    1.448    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1792(__getattr__)
       88    0.003    0.000    1.421    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1466(__init__)
  1131856    0.732    0.000    1.344    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:2453(_validate_jd_for_storage)
      464    1.311    0.003    1.311    0.003 {built-in method numpy.core._multiarray_umath._vec_string}
  3254040    0.634    0.000    1.184    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
1421322/1421250    0.467    0.000    1.163    0.000 {built-in method builtins.getattr}
   565928    0.410    0.000    1.156    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:241(jd2)
  3254040    0.621    0.000    1.064    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
     6512    0.105    0.000    1.035    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
    19952    0.951    0.000    0.994    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/indexes/base.py:5425(_getitem_slice)
30296/15276    0.034    0.000    0.952    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/decorators.py:827(__get__)
   565928    0.231    0.000    0.951    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:231(jd1)
       20    0.000    0.000    0.907    0.045 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:405(data)
       20    0.000    0.000    0.894    0.045 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:175(_get_tbdata)
13330/9062    0.009    0.000    0.792    0.000 {method 'view' of 'numpy.ndarray' objects}
       44    0.003    0.000    0.766    0.017 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:230(__array_finalize__)
5542602/5537478    0.729    0.000    0.762    0.000 {built-in method builtins.isinstance}
       44    0.019    0.000    0.712    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1524(_init_from_array)
497225/399540    0.113    0.000    0.653    0.000 {built-in method builtins.setattr}
        4    0.002    0.001    0.642    0.160 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:922(_construct_procedure)
       24    0.000    0.000    0.641    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:544(nchan)
       24    0.000    0.000    0.641    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:404(rawspectrum)
  1697832    0.443    0.000    0.620    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:441(_select_subfmts)
       24    0.000    0.000    0.596    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:834(split)
  3254056    0.551    0.000    0.551    0.000 {method 'decode' of 'bytes' objects}
   104192    0.099    0.000    0.543    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
      460    0.005    0.000    0.526    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:60(_to_bytes_or_str_array)
41084/40112    0.520    0.000    0.521    0.000 {built-in method numpy.asarray}
       20    0.000    0.000    0.508    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:396(columns)
       20    0.017    0.001    0.499    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1580(_init_from_table)
      208    0.001    0.000    0.495    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
      204    0.000    0.000    0.494    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:1162(_sorted_data)
      208    0.001    0.000    0.489    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:869(take)
      208    0.014    0.000    0.480    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:623(reindex_indexer)
   565952    0.241    0.000    0.453    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:211(in_subfmt)
  3293677    0.453    0.000    0.453    0.000 {method 'strip' of 'str' objects}
    40238    0.025    0.000    0.449    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:333(__get__)
    24812    0.031    0.000    0.444    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/card.py:282(value)
   565920    0.210    0.000    0.425    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:193(_get_allowed_subfmt)
    40238    0.099    0.000    0.424    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:442(__call__)
   565952    0.211    0.000    0.403    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:221(out_subfmt)
   282984    0.198    0.000    0.393    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:223(__len__)
     6512    0.013    0.000    0.388    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:871(name)
      944    0.112    0.000    0.385    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:255(_wrap_result)
    14080    0.018    0.000    0.351    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/blocks.py:1287(take_nd)
        4    0.000    0.000    0.347    0.087 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:314(stats)
   569906    0.339    0.000    0.339    0.000 {built-in method numpy.array}
    14284    0.018    0.000    0.318    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:59(take_nd)
    14284    0.185    0.000    0.296    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:120(_take_nd_ndarray)
       24    0.000    0.000    0.294    0.012 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:327(_str_split)
   565952    0.223    0.000    0.288    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:276(precision)
        8    0.096    0.012    0.288    0.036 /lma1/mpound/dysh/src/dysh/util/core.py:103(gbt_timestamp_to_time)
   282968    0.093    0.000    0.278    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:777(scale)
       80    0.001    0.000    0.278    0.003 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:694(__init__)
    10076    0.010    0.000    0.254    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/card.py:159(__init__)
      200    0.000    0.000    0.248    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/groupby.py:805(groups)
      200    0.000    0.000    0.248    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/ops.py:713(groups)
      200    0.001    0.000    0.248    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/groupby/grouper.py:840(groups)
   282960    0.174    0.000    0.244    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:358(<lambda>)
     4900    0.014    0.000    0.242    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/header.py:151(__getitem__)
