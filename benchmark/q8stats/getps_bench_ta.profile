using Namespace(dobench=True, key='test1', timeaverage=True, nocalibrate=False, loop=4, skipflags=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='cumulative', memory=True, quit=False)
Loading  /lma1/teuben/GBT/dysh_data/example_data/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits
STATS: {'nrows': 352, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 2, 'intnum': 11, 'sig': 1, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: getps
  name   time  VmSize VmRSS skipflags
          ms   Mbyte  Mbyte          
------- ------ ------ ----- ---------
   load  240.2 3156.2 249.0      True
getps1t 3399.6 3293.4 397.8      True
getps2t 2184.6 3232.4 337.8      True
getps3t 2043.2 3232.4 337.8      True
getps4t 2046.3 3232.4 337.8      True
 report    0.1 3232.4 337.8      True
         15458977 function calls (15003890 primitive calls) in 9.921 seconds

   Ordered by: cumulative time
   List reduced from 2858 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.001    0.000    5.401    1.350 /lma1/mpound/dysh/src/dysh/log.py:277(wrapper)
        4    0.002    0.001    5.397    1.349 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:1339(getps)
       16    0.001    0.000    4.830    0.302 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1157(__init__)
    49/13    0.001    0.000    4.488    0.345 /lma1/mpound/dysh/src/dysh/log.py:337(wrapper)
        4    0.000    0.000    4.271    1.068 /lma1/mpound/dysh/src/dysh/spectra/scan.py:638(timeaverage)
       20    0.002    0.000    3.920    0.196 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:1054(make_spectrum)
       16    0.001    0.000    3.736    0.234 /lma1/mpound/dysh/src/dysh/spectra/scan.py:505(timeaverage)
       16    0.000    0.000    3.414    0.213 /lma1/mpound/dysh/src/dysh/spectra/scan.py:191(calibrated)
       20    0.001    0.000    3.323    0.166 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:76(__init__)
       20    0.001    0.000    3.315    0.166 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/specutils/spectra/spectrum1d.py:73(__init__)
       20    0.000    0.000    3.288    0.164 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:361(pixel_to_world)
       20    0.000    0.000    3.283    0.164 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:272(values_to_high_level_objects)
       40    0.003    0.000    3.279    0.082 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:389(_get_components_and_classes)
   240/40    0.001    0.000    3.208    0.080 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/baseframe.py:1446(transform_to)
   240/40    0.002    0.000    3.207    0.080 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/composite.py:95(__call__)
   280/80    0.006    0.000    3.199    0.040 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/function.py:167(__call__)
       16    0.000    0.000    2.690    0.168 /lma1/mpound/dysh/src/dysh/spectra/scan.py:172(_finish_initialization)
       20    0.000    0.000    2.319    0.116 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:377(world_axis_object_components)
      368    0.053    0.000    2.271    0.006 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1338(exposure)
      233    0.011    0.000    2.116    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      166    0.004    0.000    2.008    0.012 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
     1568    0.005    0.000    1.859    0.001 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:124(index)
      200    0.002    0.000    1.843    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:222(itrs_to_cirs)
     5035    0.021    0.000    1.828    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4062(__getitem__)
    11620    0.183    0.000    1.790    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
      400    0.003    0.000    1.741    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/utils.py:40(get_polar_motion)
      200    0.001    0.000    1.702    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:49(cirs_to_itrs_mat)
       64    0.001    0.000    1.679    0.026 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:380(rawspectra)
86836/78923    0.053    0.000    1.669    0.000 {method 'view' of 'numpy.ndarray' objects}
     1665    0.007    0.000    1.646    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4130(_getitem_bool_array)
     1608    0.005    0.000    1.468    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6662(copy)
       83    0.005    0.000    1.466    0.018 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
     1645    0.008    0.000    1.447    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:557(copy)
       16    0.007    0.000    1.402    0.088 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1276(calibrate)
       83    0.034    0.000    1.369    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
     1000    0.001    0.000    1.357    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/state.py:52(get)
      600    0.001    0.000    1.355    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:1033(validate)
      600    0.002    0.000    1.354    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:789(open)
        1    0.000    0.000    1.325    1.325 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:610(read)
        3    0.000    0.000    1.323    0.441 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/connect.py:57(__call__)
        3    0.000    0.000    1.323    0.441 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/registry/core.py:159(read)
        3    0.000    0.000    1.306    0.435 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/connect.py:13(io_read)
        3    0.000    0.000    1.306    0.435 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/ui.py:341(read)
        2    0.000    0.000    1.298    0.649 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/cds.py:374(read)
        2    0.215    0.107    1.298    0.649 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/core.py:1382(read)
       16    0.001    0.000    1.080    0.067 /lma1/mpound/dysh/src/dysh/spectra/scan.py:489(_add_calibration_meta)
     1586    0.089    0.000    0.993    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:1782(_consolidate_inplace)
      200    0.004    0.000    0.961    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/icrs_cirs_transforms.py:71(cirs_to_icrs)
       20    0.000    0.000    0.960    0.048 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:381(world_axis_object_classes)
397474/223171    0.115    0.000    0.868    0.000 {built-in method builtins.setattr}


final 9.913971842  sec
