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
load1 12531.7 1584.7  968.9      8     66.72  533.76 16384 8064   4   2    2     20     False
load2 12382.4 1650.8 1034.2      8     66.72  533.76 16384 8064   4   2    2     20     False
load3 12296.6 1648.8 1032.6      8     66.72  533.76 16384 8064   4   2    2     20     False
load4 12287.4 1767.7 1150.8      8     66.72  533.76 16384 8064   4   2    2     20     False
         77224478 function calls (75132112 primitive calls) in 49.451 seconds

   Ordered by: cumulative time
   List reduced from 1885 to 100 due to restriction <100>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000   49.235   12.309 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.002    0.000   49.233   12.308 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
        4    0.015    0.004   44.209   11.052 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:870(_create_index_if_needed)
       32    0.037    0.001   40.463    1.264 /lma1/mpound/dysh/src/dysh/util/selection.py:1213(read)
      448    0.017    0.000   36.753    0.082 /lma1/mpound/dysh/src/dysh/util/selection.py:1062(flag)
      448    0.008    0.000   32.493    0.073 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
      448    0.005    0.000   18.267    0.041 /lma1/mpound/dysh/src/dysh/util/selection.py:446(_addrow)
      448    0.021    0.000   14.142    0.032 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
      448    0.005    0.000   13.445    0.030 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
      448    0.337    0.001   11.898    0.027 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
      448    0.053    0.000    7.782    0.017 /lma1/mpound/dysh/src/dysh/util/selection.py:419(_check_for_duplicates)
    43904    0.110    0.000    7.647    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
   281699    0.879    0.000    7.557    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:707(__array_finalize__)
      928    0.124    0.000    7.275    0.008 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3591(sort)
      452    0.004    0.000    6.862    0.015 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3169(add_row)
 1348/452    0.171    0.000    6.858    0.015 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3222(insert_row)
92480/92032    0.779    0.000    5.846    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
    44408    0.149    0.000    5.454    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1385(insert)
2918873/2018574    3.015    0.000    4.996    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1271(__setattr__)
       32    0.001    0.000    4.988    0.156 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
    28532    1.308    0.000    4.737    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:1025(indices)
    43904    0.337    0.000    4.422    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:287(get_dtypes)
    25312    0.019    0.000    4.367    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:1038(loc)
    25312    0.015    0.000    4.349    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:825(__init__)
    96392    0.445    0.000    4.299    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
   397136    4.259    0.000    4.259    0.000 {built-in method numpy.array}
   281699    0.499    0.000    4.060    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1125(_copy_attrs)
      448    0.001    0.000    4.036    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:333(_check_numbers)
42112/448    0.482    0.000    4.034    0.009 /lma1/mpound/dysh/src/dysh/util/selection.py:366(_check_type)
       32    0.003    0.000    3.942    0.123 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
1590329/1439364    0.507    0.000    3.878    0.000 {built-in method builtins.setattr}
    88608    0.299    0.000    3.434    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1363(__setitem__)
 1520/784    0.008    0.000    3.225    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
        8    0.000    0.000    3.187    0.398 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
     1488    0.004    0.000    3.173    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
464182/457586    0.246    0.000    3.163    0.000 {method 'view' of 'numpy.ndarray' objects}
        8    0.000    0.000    2.980    0.372 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
     1488    0.004    0.000    2.938    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   133280    0.274    0.000    2.863    0.000 {method 'take' of 'numpy.ndarray' objects}
    45168    0.537    0.000    2.750    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/function_base.py:5369(insert)
    25312    0.039    0.000    2.665    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:863(__getitem__)
    64520    0.022    0.000    2.540    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:252(self_iter)
    64512    0.055    0.000    2.517    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:232(__getitem__)
    25312    0.090    0.000    2.505    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:829(_get_rows)
    86752    0.205    0.000    2.496    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1339(_check_string_truncate)
    64512    0.851    0.000    2.463    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
  3172387    1.114    0.000    2.240    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:232(__get__)
  3255889    1.677    0.000    2.202    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:356(__get__)
      136    0.004    0.000    2.195    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1466(__init__)
10228403/10216863    1.759    0.000    2.135    0.000 {built-in method builtins.isinstance}
    25312    0.021    0.000    2.111    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:515(find)
    25312    0.026    0.000    2.086    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:250(find)
    25312    0.143    0.000    2.059    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/sorted_array.py:130(find)
7842146/7842074    1.427    0.000    2.046    0.000 {built-in method builtins.getattr}
   129232    1.896    0.000    1.896    0.000 {built-in method numpy.core._multiarray_umath._vec_string}
   128496    0.045    0.000    1.657    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:265(str_len)
    10064    0.163    0.000    1.592    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
    44408    0.081    0.000    1.583    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:725(__array_wrap__)
47396/23364    0.052    0.000    1.516    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/decorators.py:827(__get__)
     2408    0.008    0.000    1.478    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:503(__getitem__)
    44408    0.039    0.000    1.458    0.000 {function BaseColumn.__array_wrap__ at 0x7fa1e29e4b80}
       32    0.000    0.000    1.451    0.045 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:405(data)
       32    0.000    0.000    1.430    0.045 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:175(_get_tbdata)
      448    0.007    0.000    1.351    0.003 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/eval.py:170(eval)
    42148    0.224    0.000    1.278    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/ma/core.py:2808(__new__)
       68    0.004    0.000    1.166    0.017 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:230(__array_finalize__)
       32    0.000    0.000    1.144    0.036 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       32    0.002    0.000    1.143    0.036 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       32    0.000    0.000    1.141    0.036 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
    83328    0.243    0.000    1.129    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/ma/core.py:3217(__getitem__)
       32    0.003    0.000    1.108    0.035 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
       68    0.030    0.000    1.083    0.016 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1524(_init_from_array)
      736    0.002    0.000    1.057    0.001 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
   129056    0.213    0.000    1.000    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
    26740    0.049    0.000    0.943    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/sorted_array.py:77(_get_key_slice)
       36    0.000    0.000    0.942    0.026 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:544(nchan)
       36    0.000    0.000    0.942    0.026 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:404(rawspectrum)
       32    0.001    0.000    0.893    0.028 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:75(_init_flags)
     1508    0.247    0.000    0.892    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
134749/131389    0.090    0.000    0.869    0.000 {built-in method builtins.all}
    83812    0.160    0.000    0.855    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/ma/core.py:2978(__array_finalize__)
   161024    0.153    0.000    0.834    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
   205137    0.190    0.000    0.816    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:582(__set__)
       32    0.000    0.000    0.810    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/hdu/table.py:396(columns)
       32    0.025    0.001    0.795    0.025 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:1580(_init_from_table)
    24864    0.035    0.000    0.783    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:1445(equals)
      896    0.001    0.000    0.775    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:540(insert_row)
      896    0.081    0.000    0.773    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:164(insert_row)
    91992    0.219    0.000    0.769    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:1863(from_array)
      928    0.002    0.000    0.764    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:730(__exit__)
2915531/2433458    0.554    0.000    0.757    0.000 {built-in method builtins.len}
   125476    0.443    0.000    0.755    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/ma/core.py:2952(_update_from)
      508    0.003    0.000    0.741    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
    24864    0.075    0.000    0.741    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/base.py:144(equals)
    64974    0.041    0.000    0.731    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:333(__get__)
   104272    0.194    0.000    0.727    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:6301(__setattr__)
       32    0.000    0.000    0.725    0.023 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:501(drop_fields)
       32    0.179    0.006    0.723    0.023 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:35(recursive_fill_fields)
      448    0.001    0.000    0.705    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/expr.py:796(__init__)
      448    0.001    0.000    0.703    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/expr.py:824(parse)
