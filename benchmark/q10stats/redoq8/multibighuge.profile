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
load1 36955.3 3761.3 2050.8      3    768.09 2304.27 32768 17496   3   2    2    102     False
load2 36838.4 3876.7 2165.9      3    768.09 2304.27 32768 17496   3   2    2    102     False
load3 36835.6 3878.1 2167.7      3    768.09 2304.27 32768 17496   3   2    2    102     False
load4 37031.2 3890.7 2179.6      3    768.09 2304.27 32768 17496   3   2    2    102     False
         241156520 function calls (234961554 primitive calls) in 147.520 seconds

   Ordered by: cumulative time
   List reduced from 1901 to 100 due to restriction <100>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.000    0.000  147.225   36.806 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.001    0.000  147.222   36.806 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:75(__init__)
        4    0.016    0.004  142.439   35.610 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:870(_create_index_if_needed)
       12    0.039    0.003  134.357   11.196 /lma1/mpound/dysh/src/dysh/util/selection.py:1213(read)
     1152    0.023    0.000  121.931    0.106 /lma1/mpound/dysh/src/dysh/util/selection.py:1062(flag)
     1152    0.022    0.000  119.871    0.104 /lma1/mpound/dysh/src/dysh/util/selection.py:489(_base_select)
     1152    0.012    0.000   82.204    0.071 /lma1/mpound/dysh/src/dysh/util/selection.py:446(_addrow)
     1152    0.350    0.000   51.567    0.045 /lma1/mpound/dysh/src/dysh/util/selection.py:419(_check_for_duplicates)
     1152    0.054    0.000   37.445    0.033 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4680(query)
     1152    0.013    0.000   35.604    0.031 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4846(eval)
     1152    0.864    0.001   30.856    0.027 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:644(_get_cleaned_column_resolvers)
   174572    8.171    0.000   29.563    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:1025(indices)
   166464    0.119    0.000   28.758    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:1038(loc)
   166464    0.095    0.000   28.639    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:825(__init__)
     2316    0.317    0.000   25.125    0.011 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3591(sort)
   918519    2.891    0.000   24.780    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:707(__array_finalize__)
   112928    0.430    0.000   19.905    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:6432(dtypes)
     1156    0.010    0.000   17.818    0.015 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3169(add_row)
3460/1156    0.445    0.000   17.807    0.015 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/table.py:3222(insert_row)
   166464    0.260    0.000   17.532    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:863(__getitem__)
   166464    0.596    0.000   16.490    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:829(_get_rows)
9329921/6436994    9.694    0.000   16.138    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1271(__setattr__)
   221172    0.812    0.000   15.356    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1363(__setitem__)
234588/233436    1.991    0.000   14.933    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/series.py:389(__init__)
   114104    0.389    0.000   14.144    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1385(insert)
   166464    0.139    0.000   13.886    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:515(find)
   166464    0.168    0.000   13.717    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:250(find)
   918519    1.661    0.000   13.612    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1125(_copy_attrs)
   166464    0.927    0.000   13.549    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/sorted_array.py:130(find)
 18663463    6.551    0.000   13.052    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:232(__get__)
   216540    0.527    0.000   11.802    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1339(_check_string_truncate)
   112928    0.917    0.000   11.421    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:287(get_dtypes)
 18172537    8.418    0.000   11.225    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:356(__get__)
   683750   10.653    0.000   10.653    0.000 {built-in method numpy.array}
4166121/4103956    1.426    0.000   10.355    0.000 {built-in method builtins.setattr}
   324736   10.181    0.000   10.181    0.000 {built-in method numpy.core._multiarray_umath._vec_string}
   324460    0.113    0.000    9.664    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/core/defchararray.py:265(str_len)
   236236    1.095    0.000    9.585    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/construction.py:517(sanitize_array)
   339376    0.828    0.000    7.294    0.000 {method 'take' of 'numpy.ndarray' objects}
   114424    1.431    0.000    7.152    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/numpy/lib/function_base.py:5369(insert)
28854446/28854374    4.987    0.000    6.734    0.000 {built-in method builtins.getattr}
        8    0.000    0.000    6.400    0.800 /lma1/mpound/dysh/src/dysh/util/selection.py:55(__init__)
        8    0.000    0.000    6.302    0.788 /lma1/mpound/dysh/src/dysh/util/selection.py:94(_add_utc_column)
  600/324    0.003    0.000    6.300    0.019 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4271(__setitem__)
      588    0.002    0.000    6.279    0.011 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:4514(_set_item)
      588    0.002    0.000    6.174    0.010 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:5242(_sanitize_column)
   170356    0.305    0.000    5.977    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/sorted_array.py:77(_get_key_slice)
26045568/26030624    4.540    0.000    5.371    0.000 {built-in method builtins.isinstance}
   139976    0.048    0.000    5.353    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:252(self_iter)
   139968    0.121    0.000    5.305    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/shapes.py:232(__getitem__)
   139968    1.719    0.000    5.184    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
   165312    0.249    0.000    4.975    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:1445(equals)
424725/416597    0.279    0.000    4.946    0.000 {built-in method builtins.all}
       12    0.000    0.000    4.690    0.391 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:37(__init__)
   165312    0.494    0.000    4.684    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/base.py:144(equals)
     1152    0.018    0.000    4.217    0.004 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/eval.py:170(eval)
       12    0.006    0.000    4.083    0.340 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:150(create_index)
   114104    0.211    0.000    4.059    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:725(__array_wrap__)
   495936    0.259    0.000    4.028    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/base.py:155(<genexpr>)
 18887071    3.858    0.000    3.858    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:335(_parent)
   332184    0.687    0.000    3.776    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/indexes/base.py:5552(equals)
   114104    0.092    0.000    3.738    0.000 {function BaseColumn.__array_wrap__ at 0x7fabc8b40b80}
889930/887214    0.445    0.000    3.511    0.000 {method 'view' of 'numpy.ndarray' objects}
   338408    0.247    0.000    2.907    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/sorted_array.py:5(_searchsorted)
 18635414    2.902    0.000    2.902    0.000 {method 'get' of 'dict' objects}
     1152    0.005    0.000    2.796    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/engines.py:65(evaluate)
   621685    0.603    0.000    2.588    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:582(__set__)
   338408    0.728    0.000    2.578    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1001(searchsorted)
   221172    0.873    0.000    2.255    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:617(adjust_indices)
8269275/6771258    1.612    0.000    2.226    0.000 {built-in method builtins.len}
   279968    0.458    0.000    2.145    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
   918519    0.758    0.000    2.038    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:767(name)
     2304    0.004    0.000    2.030    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:540(insert_row)
     2304    0.209    0.000    2.025    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:164(insert_row)
   233376    0.552    0.000    1.995    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:1863(from_array)
   621685    1.279    0.000    1.984    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:371(__set__)
     1348    0.009    0.000    1.966    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:4027(take)
     2316    0.005    0.000    1.926    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/index.py:730(__exit__)
      568    0.483    0.001    1.908    0.003 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
       12    0.000    0.000    1.842    0.153 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/frame.py:10170(apply)
       12    0.003    0.000    1.841    0.153 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:864(apply)
       12    0.000    0.000    1.838    0.153 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1061(apply_standard)
       12    0.001    0.000    1.819    0.152 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/apply.py:1070(apply_series_generator)
      276    0.001    0.000    1.796    0.007 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:173(<lambda>)
     1184    0.004    0.000    1.785    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/indexing.py:1176(__getitem__)
     1152    0.004    0.000    1.777    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/indexing.py:1397(_getitem_axis)
   247176    0.462    0.000    1.768    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:6301(__setattr__)
   152088    0.444    0.000    1.760    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/missing.py:466(array_equivalent)
     1152    0.003    0.000    1.753    0.002 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/indexing.py:1205(_getbool_axis)
   117504    0.133    0.000    1.726    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/parsing.py:99(clean_column_name)
   227632    0.149    0.000    1.674    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:511(_validate_dtype)
     1152    0.009    0.000    1.672    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/generic.py:4142(_take_with_is_copy)
   114104    0.314    0.000    1.658    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:129(_expand_string_array_for_values)
     1348    0.006    0.000    1.647    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:869(take)
     1348    0.061    0.000    1.574    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:623(reindex_indexer)
   324460    0.624    0.000    1.559    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/utils/misc.py:554(dtype_bytes_or_chars)
     1152    0.005    0.000    1.552    0.001 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/computation/engines.py:112(_evaluate)
   519404    1.540    0.000    1.540    0.000 {method 'reduce' of 'numpy.ufunc' objects}
  1035741    0.471    0.000    1.534    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:599(parent_table)
   228856    0.638    0.000    1.528    0.000 /n/lma1/mpound/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/common.py:1596(pandas_dtype)


