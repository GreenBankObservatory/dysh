using Namespace(dobench=True, key='test1', timeaverage=False, nocalibrate=False, loop=4, feeds=16, skipflags=True, big=False, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='time', memory=True, quit=False)
Loading  /lma1/teuben/GBT/dysh_data/example_data/mapping-L/data/TGBT17A_506_11.raw.vegas
STATS: {'nrows': 15680, 'nfiles': 1, 'fdnum': 1, 'ifnum': 5, 'plnum': 2, 'intnum': 61, 'sig': 2, 'cal': 2, 'nchan': 32768}
Found ostype= linux
Removing file.fits from ngc6946
 ID    TAG                  SCAN               IFNUM PLNUM FDNUM # SELECTED
--- --------- -------------------------------- ----- ----- ----- ----------
  0 f7e78b2f6 [14,15,16,17,18...3,24,25,26,27]     0     0     0       1708
STATS: {'nrows': 1708, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 1, 'intnum': 61, 'sig': 1, 'cal': 2, 'nchan': 4096}
# Dysh Benchmark: otf
   name     time  VmSize VmRSS  skipflags
             ms   Mbyte  Mbyte           
--------- ------- ------ ------ ---------
    load1  5222.7 3860.5 2260.1      True
   write1  2541.4 4949.8 3350.2      True
    load2   556.6 4984.2 3384.7      True
getsigref 18536.5 1444.6  828.7      True
   write2   949.8 1492.5  876.7      True
   report     0.1 1492.5  876.7      True
         47616860 function calls (45189442 primitive calls) in 27.772 seconds

   Ordered by: internal time
   List reduced from 3245 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
1854872/15516    1.671    0.000    5.588    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:118(deepcopy)
   147965    1.586    0.000    5.316    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
5330553/5325435    0.797    0.000    0.892    0.000 {built-in method builtins.isinstance}
    99613    0.695    0.000    0.698    0.000 {method 'copy' of 'numpy.ndarray' objects}
  4272896    0.657    0.000    0.657    0.000 {method 'get' of 'dict' objects}
609578/603632    0.626    0.000    1.101    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1792(__getattr__)
    11859    0.543    0.000    0.988    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:2276(_merge_blocks)
   301092    0.506    0.000    0.539    0.000 {built-in method numpy.array}
   298693    0.499    0.000    2.434    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
2035211/1938589    0.497    0.000    1.442    0.000 {built-in method builtins.getattr}
     5825    0.467    0.000    0.467    0.000 {built-in method numpy.core._multiarray_umath._vec_string}
   597580    0.412    0.000    0.717    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:2453(_validate_jd_for_storage)
     4352    0.409    0.000    4.648    0.001 {method '__deepcopy__' of 'numpy.ndarray' objects}
     6786    0.298    0.000    0.358    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/numpy/core/shape_base.py:219(vstack)
303955/211896    0.279    0.000    0.518    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1271(__setattr__)
  1812272    0.266    0.000    0.266    0.000 {method 'strip' of 'str' objects}
     1784    0.253    0.000    1.801    0.001 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/internals/managers.py:1782(_consolidate_inplace)
    42084    0.252    0.000    0.578    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/ascii/fixedwidth.py:36(__call__)
   895366    0.237    0.000    0.334    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:441(_select_subfmts)
      142    0.234    0.002    0.987    0.007 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:46(_str_map)
41734/39726    0.228    0.000    0.276    0.000 {built-in method numpy.asarray}
  1470791    0.228    0.000    0.228    0.000 {method 'append' of 'list' objects}
  1988034    0.219    0.000    0.219    0.000 {built-in method builtins.id}
   298790    0.218    0.000    0.614    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:241(jd2)
        2    0.213    0.106    1.225    0.612 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/ascii/core.py:1380(read)
      780    0.211    0.000    3.529    0.005 /lma1/teuben/GBT/dysh/src/dysh/spectra/scan.py:1335(exposure)
   172450    0.209    0.000    0.209    0.000 {built-in method builtins.locals}
     2136    0.205    0.000    0.209    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/indexes/base.py:5425(_getitem_slice)
        3    0.202    0.067    0.649    0.216 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/numpy/lib/recfunctions.py:35(recursive_fill_fields)
     7271    0.191    0.000    0.195    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/fits/header.py:1811(_updateindices)
   941406    0.189    0.000    0.326    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/ascii/core.py:405(process_val)
    11488    0.186    0.000    1.697    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
   183808    0.172    0.000    0.776    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
    73164    0.171    0.000    0.190    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/utils/data_info.py:356(__get__)
614977/442654    0.166    0.000    1.106    0.000 {built-in method builtins.setattr}
55977/55899    0.161    0.000    0.162    0.000 {method 'reduce' of 'numpy.ufunc' objects}
    59068    0.158    0.000    0.658    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/config/configuration.py:442(__call__)
   792764    0.154    0.000    0.295    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/strings/accessor.py:2002(<lambda>)
    18368    0.151    0.000    0.327    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:120(_take_nd_ndarray)
   792764    0.150    0.000    0.261    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:451(<lambda>)
      117    0.144    0.001    0.538    0.005 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:1580(construct_1d_object_array_from_listlike)
875294/807586    0.143    0.000    0.171    0.000 {built-in method builtins.len}
   795632    0.141    0.000    0.141    0.000 {method 'decode' of 'bytes' objects}
  1776126    0.138    0.000    0.138    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:172(_deepcopy_atomic)
        4    0.138    0.034    0.500    0.125 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:281(from_columns)
    38407    0.129    0.000    0.223    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/numpy/ma/core.py:2952(_update_from)
   298693    0.129    0.000    0.245    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:211(in_subfmt)
   298693    0.124    0.000    0.160    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:276(precision)
    68936    0.124    0.000    0.140    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/strings/object_array.py:358(<lambda>)
   298693    0.124    0.000    0.227    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:221(out_subfmt)


final 27.806951273  sec
