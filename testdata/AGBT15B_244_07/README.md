# AGBT15B_244_07

Test data from project AGBT15B_244

The test data contains
[Nod](https://gbtdocs.readthedocs.io/en/latest/references/software/scheduling-blocks.html#astrid_commands.Nod)
observations of M82.

## Data preparation

The original data is too large (729MB) for testing, taken in 4 IFs and
2 PL's.   Two calibration scans were done at the start and end of
the science observations (5 nods). The calibration scan was 41 integrations, the Nod
consisted of 61 integrations.

The small test dataset is just 3 integrations (in order to capture the
sky, cold1 and cold2 settings of the calibration wheel), and retrieved
as follows from the large dataset. See test_vane.py for runnable code.

```
  sdf3=GBTFITSLoad(dysh_data(accept='AGBT15B_244_07/AGBT15B_244_07.raw.vegas'))
  sdf3.summary()

  mkdir("AGBT15B_244_07_test")
  intnums=[1,14,31]
  sdf3.write("AGBT15B_244_07_test/file.fits", scan=[130,131,132], ifnum=1, plnum=0,  intnum=intnums, overwrite=True)

  test3 = GBTFITSLoad("AGBT15B_244_07_test")    
  test3.summary()


   SCAN OBJECT VELOCITY    PROC  PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0   130    M82      0.0  CALSEQ         1    87.23    86.4    1     1     3      2  334.378168  46.559458
1   131    M82      0.0     Nod         1    87.23    86.4    1     1     3      2  334.347874  46.476138
2   132    M82      0.0     Nod         2    87.23    86.4    1     1     3      2  334.427945  46.377163


```

No tests comparing with GBTIDL were done yet.
