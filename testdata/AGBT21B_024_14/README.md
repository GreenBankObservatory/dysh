# GBT21B_024_14

Test data from project GBT21B-024

The test data contains [Nod](https://gbtdocs.readthedocs.io/en/latest/references/software/scheduling-blocks.html#astrid_commands.Nod) observations of [NGC5908](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC5908),
as a test in the EDGE survey to check the flux between previous CARMA obvious and the new ARGUS observations. In addition a VANE/SKY calibration scan was done, from which the Tsys can be determined to scale the
Nod observations.

Tcal was determined to be 272 K using the vanecal.pro online. In case dysh is being run offsite, this value can be entered manually.


## Data preparation
The original data is large (1640MB), and consists of almost six hours of observing. Towards the end NGC5908 was observed in a Nod procedure, preceded by a VANE/SKY calibration. The calibration was 21 integrations,
the Nod consisted of 61 integrations.

The small test dataset is just 1 integration, and retrieved as follows from the large dataset.

```
   sdf2=GBTFITSLoad(dysh_data('AGBT21B_024_14/AGBT21B_024_14.raw.vegas'), skipflags=True)
   sdf2.summary()

   mkdir("AGBT21B_024_14_test")
   sdf2.write("AGBT21B_024_14_test/file.fits", scan=range(329,335), intnum=0, overwrite=True)

   test2 = GBTFITSLoad("AGBT21B_024_14_test")
   test2.summary()
   
   SCAN   OBJECT VELOCITY   PROC  PROCSEQN   RESTFREQ    DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0   329     VANE      0.0  Track         1  114.03043  114.03043    1     1     1     16  333.863964  70.214864
1   330      SKY      0.0  Track         1  114.03043  114.03043    1     1     1     16  333.856661  70.214412
2   331  NGC5908      0.0    Nod         1  114.03043  114.03043    1     1     1     16  333.171573  70.123804
3   332  NGC5908      0.0    Nod         2  114.03043  114.03043    1     1     1     16  332.950559   70.04885
4   333  NGC5908      0.0    Nod         1  114.03043  114.03043    1     1     1     16  332.734249  70.008266
5   334  NGC5908      0.0    Nod         2  114.03043  114.03043    1     1     1     16  332.519705  69.928914
```

## Bugs

1) The test data does not cleanly show a summary in gbtidl, instead shows each file.
2) No index, it was created
3) vanecal.pro runs into infinite loop
