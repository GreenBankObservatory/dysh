using Namespace(dobench=True, key='test1', timeaverage=True, nocalibrate=False, loop=4, skipflags=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='cumulative', memory=True, quit=False)
Loading  /lma1/teuben/GBT/dysh_data/example_data/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits
STATS: {'nrows': 352, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 2, 'intnum': 11, 'sig': 1, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: getps
  name   time  VmSize VmRSS skipflags
          ms   Mbyte  Mbyte          
------- ------ ------ ----- ---------
   load  241.1 3156.1 248.8      True
getps1t 3417.8 3293.5 398.1      True
getps2t 2185.0 3231.3 336.9      True
getps3t 2054.0 3231.3 336.9      True
getps4t 2044.2 3231.3 336.9      True
 report    0.1 3231.3 336.9      True
         15459296 function calls (15004209 primitive calls) in 9.949 seconds

   Ordered by: cumulative time
   List reduced from 2858 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        4    0.001    0.000    5.418    1.354 /lma1/mpound/dysh/src/dysh/log.py:278(wrapper)
        4    0.002    0.001    5.415    1.354 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:1339(getps)
       16    0.001    0.000    4.850    0.303 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1154(__init__)
    49/13    0.001    0.000    4.499    0.346 /lma1/mpound/dysh/src/dysh/log.py:338(wrapper)
        4    0.000    0.000    4.281    1.070 /lma1/mpound/dysh/src/dysh/spectra/scan.py:638(timeaverage)
       20    0.002    0.000    3.935    0.197 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:1054(make_spectrum)
       16    0.001    0.000    3.751    0.234 /lma1/mpound/dysh/src/dysh/spectra/scan.py:505(timeaverage)
       16    0.000    0.000    3.432    0.214 /lma1/mpound/dysh/src/dysh/spectra/scan.py:191(calibrated)
       20    0.001    0.000    3.342    0.167 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:76(__init__)
       20    0.001    0.000    3.334    0.167 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/specutils/spectra/spectrum1d.py:73(__init__)
       20    0.000    0.000    3.308    0.165 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:361(pixel_to_world)
       20    0.000    0.000    3.303    0.165 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:272(values_to_high_level_objects)
       40    0.003    0.000    3.299    0.082 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:389(_get_components_and_classes)
   240/40    0.001    0.000    3.228    0.081 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/baseframe.py:1446(transform_to)
   240/40    0.002    0.000    3.226    0.081 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/composite.py:95(__call__)
   280/80    0.006    0.000    3.218    0.040 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/function.py:167(__call__)
       16    0.000    0.000    2.694    0.168 /lma1/mpound/dysh/src/dysh/spectra/scan.py:172(_finish_initialization)
       20    0.000    0.000    2.335    0.117 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:377(world_axis_object_components)
      368    0.053    0.000    2.273    0.006 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1335(exposure)
      233    0.011    0.000    2.128    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      166    0.004    0.000    2.019    0.012 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
     1568    0.005    0.000    1.857    0.001 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:124(index)
      200    0.002    0.000    1.853    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:222(itrs_to_cirs)
     5035    0.021    0.000    1.827    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4062(__getitem__)
    11620    0.189    0.000    1.800    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
      400    0.003    0.000    1.757    0.004 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/utils.py:40(get_polar_motion)
      200    0.001    0.000    1.710    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:49(cirs_to_itrs_mat)
       64    0.001    0.000    1.687    0.026 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:380(rawspectra)
86836/78923    0.052    0.000    1.666    0.000 {method 'view' of 'numpy.ndarray' objects}
     1665    0.007    0.000    1.643    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/frame.py:4130(_getitem_bool_array)
     1608    0.005    0.000    1.469    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/generic.py:6662(copy)
       83    0.005    0.000    1.462    0.018 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
     1645    0.008    0.000    1.447    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:557(copy)
       16    0.007    0.000    1.406    0.088 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1273(calibrate)
     1000    0.001    0.000    1.370    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/state.py:52(get)
      600    0.001    0.000    1.368    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:1033(validate)
      600    0.002    0.000    1.368    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:789(open)
       83    0.034    0.000    1.365    0.016 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
        1    0.000    0.000    1.338    1.338 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:610(read)
        3    0.000    0.000    1.337    0.446 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/connect.py:57(__call__)
        3    0.000    0.000    1.337    0.446 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/registry/core.py:159(read)
        3    0.000    0.000    1.320    0.440 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/connect.py:13(io_read)
        3    0.000    0.000    1.320    0.440 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/ui.py:341(read)
        2    0.000    0.000    1.312    0.656 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/cds.py:374(read)
        2    0.216    0.108    1.312    0.656 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/core.py:1382(read)
       16    0.001    0.000    1.080    0.068 /lma1/mpound/dysh/src/dysh/spectra/scan.py:489(_add_calibration_meta)
     1586    0.088    0.000    0.991    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/pandas/core/internals/managers.py:1782(_consolidate_inplace)
       20    0.000    0.000    0.965    0.048 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:381(world_axis_object_classes)
      200    0.004    0.000    0.964    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/icrs_cirs_transforms.py:71(cirs_to_icrs)
397474/223171    0.116    0.000    0.876    0.000 {built-in method builtins.setattr}


final 9.942146704  sec
