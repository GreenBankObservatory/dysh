using Namespace(dobench=True, key='test1', timeaverage=False, nocalibrate=False, loop=4, feeds=16, skipflags=True, big=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='time', memory=True, quit=False)
Loading  /lma1/teuben/GBT/dysh_data/example_data/mapping-L/data/TGBT17A_506_11.raw.vegas
STATS: {'nrows': 15680, 'nfiles': 1, 'fdnum': 1, 'ifnum': 5, 'plnum': 2, 'intnum': 61, 'sig': 2, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: otf
   name     time   VmSize VmRSS  skipflags
             ms    Mbyte  Mbyte           
--------- -------- ------ ------ ---------
    load1   5241.1 3860.7 2260.0      True
   write1      0.1 3860.7 2260.0      True
    load2      0.0 3860.7 2260.0      True
getsigref 115682.9 3947.3 2348.3      True
   write2    930.9 3950.7 2351.8      True
   report      0.1 3950.7 2351.8      True
         266045053 function calls (233077953 primitive calls) in 121.854 seconds

   Ordered by: internal time
   List reduced from 3188 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
32442609/15106   28.287    0.000   86.355    0.006 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:118(deepcopy)
   919141   11.354    0.000   35.758    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1350(_apply)
   231338   10.629    0.000   13.679    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:120(_take_nd_ndarray)
 66432652    9.800    0.000    9.800    0.000 {method 'get' of 'dict' objects}
     4589    6.348    0.001   85.474    0.019 {method '__deepcopy__' of 'numpy.ndarray' objects}
 33382677    3.685    0.000    3.685    0.000 {built-in method builtins.id}
  1841032    3.118    0.000   14.663    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:155(__init__)
21853387/21848613    2.984    0.000    3.093    0.000 {built-in method builtins.isinstance}
  3682276    2.941    0.000    4.518    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:2453(_validate_jd_for_storage)
      780    2.776    0.004   18.350    0.024 /lma1/teuben/GBT/dysh/src/dysh/spectra/scan.py:1335(exposure)
 31556088    2.422    0.000    2.422    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:172(_deepcopy_atomic)
2043383/2037437    2.387    0.000    2.871    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1792(__getattr__)
5779720/5691410    1.640    0.000    4.447    0.000 {built-in method builtins.getattr}
  5522393    1.494    0.000    2.110    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:441(_select_subfmts)
  1841138    1.385    0.000    3.902    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:241(jd2)
   236396    1.298    0.000    1.298    0.000 {built-in method numpy.empty}
  1841032    0.890    0.000    1.128    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:276(precision)
  1841032    0.812    0.000    1.537    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:211(in_subfmt)
  1841138    0.771    0.000    3.184    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:231(jd1)
  1841032    0.725    0.000    1.377    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:221(out_subfmt)
  1839063    0.713    0.000    1.445    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:193(_get_allowed_subfmt)
   851799    0.674    0.000   33.997    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1294(copy)
   229674    0.659    0.000    0.708    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/dtypes/cast.py:551(maybe_promote)
   852674    0.650    0.000    0.934    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:61(copy)
   882613    0.517    0.000    0.836    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/copy.py:231(_keep_alive)
   227437    0.495    0.000    0.528    0.000 {built-in method numpy.array}
  3770308    0.487    0.000    0.487    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:227(jd1)
   924054    0.437    0.000    0.677    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1097(masked)
   231338    0.428    0.000    1.324    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:564(_take_preprocess_indexer_and_fill_value)
  1841451    0.412    0.000    0.412    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:2469(_broadcast_writeable)
   851799    0.405    0.000   34.402    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:1472(__deepcopy__)
   931355    0.384    0.000    0.616    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:262(scale)
  2085559    0.357    0.000    0.357    0.000 {method 'append' of 'list' objects}
   231263    0.345    0.000   14.713    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/internals/blocks.py:1287(take_nd)
     5560    0.333    0.000    0.333    0.000 {built-in method numpy.core._multiarray_umath._vec_string}
   930424    0.318    0.000    0.933    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:777(scale)
   231338    0.317    0.000   14.067    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:59(take_nd)
98538/98460    0.317    0.000    0.318    0.000 {method 'reduce' of 'numpy.ufunc' objects}
1360012/1198804    0.312    0.000    1.196    0.000 {built-in method builtins.setattr}
   920441    0.299    0.000    0.419    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:845(precision)
   922133    0.296    0.000    0.415    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:871(out_subfmt)
   920441    0.287    0.000    0.399    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/core.py:858(in_subfmt)
  1841032    0.285    0.000    0.285    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:268(scale)
292995/204257    0.270    0.000    0.499    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/table/column.py:1271(__setattr__)
     6951    0.265    0.000    0.288    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/array_algos/take.py:353(wrapper)
    42084    0.259    0.000    0.587    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/io/ascii/fixedwidth.py:36(__call__)
  1767160    0.258    0.000    0.258    0.000 {method 'strip' of 'str' objects}
   231308    0.252    0.000    0.252    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/pandas/core/internals/blocks.py:292(make_block_same_class)
  1843081    0.252    0.000    0.252    0.000 /lma1/teuben/GBT/anaconda3/lib/python3.12/site-packages/astropy/time/formats.py:237(jd2)
   974279    0.252    0.000    0.252    0.000 {built-in method __new__ of type object at 0x936bc0}


final 121.855038186  sec
