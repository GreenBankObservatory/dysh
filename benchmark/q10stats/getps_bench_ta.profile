using Namespace(dobench=True, key='test1', timeaverage=True, nocalibrate=False, loop=4, skipflags=True, out=None, append=False, overwrite=False, profile=True, statslines='50', sortkey='cumulative', memory=True, quit=False)
AGBT05B_047_01.raw.acs.fits already downloaded
Loading  AGBT05B_047_01.raw.acs.fits
STATS: {'nrows': 352, 'nfiles': 1, 'fdnum': 1, 'ifnum': 1, 'plnum': 2, 'intnum': 11, 'sig': 1, 'cal': 2, 'nchan': 32768}
Found ostype= linux
# Dysh Benchmark: getps
  name   time  VmSize VmRSS skipflags
          ms   Mbyte  Mbyte          
------- ------ ------ ----- ---------
   load  249.5 3156.2 249.5      True
getps1t 2969.2 3263.4 368.2      True
getps2t 1649.4 3230.4 335.6      True
getps3t 1514.7 3230.4 335.6      True
getps4t 1509.9 3230.4 335.6      True
 report    0.1 3230.4 335.6      True
         12670600 function calls (12242494 primitive calls) in 7.900 seconds

   Ordered by: cumulative time
   List reduced from 2924 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    49/13    0.001    0.000    4.373    0.336 /lma1/mpound/dysh/src/dysh/log.py:342(wrapper)
        4    0.000    0.000    4.149    1.037 /lma1/mpound/dysh/src/dysh/spectra/scan.py:950(timeaverage)
       20    0.002    0.000    4.009    0.200 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:1106(make_spectrum)
       16    0.001    0.000    3.619    0.226 /lma1/mpound/dysh/src/dysh/spectra/scan.py:761(timeaverage)
       16    0.000    0.000    3.508    0.219 /lma1/mpound/dysh/src/dysh/spectra/scan.py:335(calibrated)
        4    0.001    0.000    3.492    0.873 /lma1/mpound/dysh/src/dysh/log.py:282(wrapper)
        4    0.003    0.001    3.488    0.872 /lma1/mpound/dysh/src/dysh/fits/gbtfitsload.py:1388(getps)
       20    0.001    0.000    3.399    0.170 /lma1/mpound/dysh/src/dysh/spectra/spectrum.py:80(__init__)
       20    0.001    0.000    3.391    0.170 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/specutils/spectra/spectrum1d.py:73(__init__)
       20    0.000    0.000    3.368    0.168 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:361(pixel_to_world)
       20    0.000    0.000    3.363    0.168 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/high_level_api.py:272(values_to_high_level_objects)
       40    0.003    0.000    3.359    0.084 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:389(_get_components_and_classes)
   240/40    0.001    0.000    3.288    0.082 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/baseframe.py:1446(transform_to)
   240/40    0.002    0.000    3.286    0.082 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/composite.py:95(__call__)
   280/80    0.006    0.000    3.278    0.041 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/transformations/function.py:167(__call__)
       16    0.001    0.000    2.889    0.181 /lma1/mpound/dysh/src/dysh/spectra/scan.py:1499(__init__)
       20    0.000    0.000    2.408    0.120 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:377(world_axis_object_components)
      233    0.011    0.000    2.144    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:506(__getitem__)
      166    0.005    0.000    2.036    0.012 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1467(__init__)
      200    0.002    0.000    1.907    0.010 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:222(itrs_to_cirs)
      400    0.003    0.000    1.831    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/utils.py:40(get_polar_motion)
      200    0.001    0.000    1.760    0.009 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/intermediate_rotation_transforms.py:49(cirs_to_itrs_mat)
       64    0.001    0.000    1.708    0.027 /lma1/mpound/dysh/src/dysh/fits/sdfitsload.py:380(rawspectra)
84036/76143    0.049    0.000    1.682    0.000 {method 'view' of 'numpy.ndarray' objects}
    11620    0.185    0.000    1.681    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:561(__init__)
       83    0.005    0.000    1.490    0.018 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/fitsrec.py:231(__array_finalize__)
     1000    0.001    0.000    1.420    0.001 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/state.py:52(get)
      600    0.001    0.000    1.419    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:1033(validate)
      600    0.002    0.000    1.418    0.002 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:789(open)
        1    0.000    0.000    1.388    1.388 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:610(read)
        3    0.000    0.000    1.386    0.462 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/table/connect.py:57(__call__)
        3    0.000    0.000    1.386    0.462 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/registry/core.py:159(read)
        3    0.000    0.000    1.369    0.456 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/connect.py:13(io_read)
        3    0.000    0.000    1.368    0.456 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/ui.py:341(read)
        2    0.000    0.000    1.360    0.680 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/cds.py:374(read)
        2    0.215    0.107    1.360    0.680 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/core.py:1382(read)
       83    0.034    0.000    1.258    0.015 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1525(_init_from_array)
      200    0.004    0.000    0.983    0.005 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/builtin_frames/icrs_cirs_transforms.py:71(cirs_to_icrs)
       20    0.000    0.000    0.951    0.048 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/wcs/wcsapi/fitswcs.py:381(world_axis_object_classes)
397474/223171    0.113    0.000    0.872    0.000 {built-in method builtins.setattr}
   185920    0.175    0.000    0.754    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:532(__set__)
       16    0.001    0.000    0.722    0.045 /lma1/mpound/dysh/src/dysh/spectra/scan.py:301(_finish_initialization)
        1    0.000    0.000    0.629    0.629 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:559(_combine_a_b_columns)
        1    0.000    0.000    0.625    0.625 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:964(_substitute_iers_b)
      200    0.003    0.000    0.623    0.003 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/coordinates/erfa_astrom.py:38(apco)
        1    0.000    0.000    0.623    0.623 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:230(open)
        1    0.000    0.000    0.623    0.623 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/utils/iers/iers.py:708(read)
666540/609014    0.163    0.000    0.610    0.000 {built-in method builtins.getattr}
    42524    0.266    0.000    0.596    0.000 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/ascii/fixedwidth.py:36(__call__)
       82    0.004    0.000    0.564    0.007 /lma1/mpound/.local/share/hatch/env/virtual/dysh/rN9ljKHD/dysh/lib/python3.12/site-packages/astropy/io/fits/column.py:1512(_init_from_coldefs)


final 7.89285912  sec
