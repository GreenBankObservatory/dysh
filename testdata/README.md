### Test data directory

Here is a summary, but each directory has its own explanation in detail.

20m/
AGBT04A_008_02/             testtrim; multiple binary tables
AGBT05B_047_01/             getps   NGC5291 11 int, 1420 MHz
AGBT13A_124_06/
AGBT15B_244_07/
AGBT17A_404_01/                     A123606
AGBT17B_173_04/
AGBT17B_456_03/
AGBT18A_333_21/
AGBT18B_354_03/
AGBT20B_014_03.raw.vegas/   getfs   M33S    [scan*34][int*11][if*4][pol*2][sig*2][cal*2]  has 11840/11968 rows
AGBT20B_295_02
AGBT21B_024_01
AGBT21B_024_14
AGBT22A_325_15
calibration
extract_rows.py
gbtidl_spectra
gshift_box.fits
TGBT17A_506_11
TGBT21A_501_11/             getps   NGC2415
TGBT21A_501_11/             gettp
               TGBT21A_501_11.raw.vegas.fits : Nrows: 6040   Ncols: 74  Nchan: 32768  nP: 791.675 Mp
TGBT21A_504_01/             getfs   W3OH    11 int, 2 pol, 1660 MHz.
TGBT22A_503_02              getfs   W3_1       test=    example=
TRCO_230413_Ka              subbeamnod
TSCAL_220105_W


----

old/removed:

TREG_050627.raw.acs.fits    getfs   W3OH,IC1481  

TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits   File with two binary tables of length 3 row and 5 rows.  Good to test multi-bintable functions.
