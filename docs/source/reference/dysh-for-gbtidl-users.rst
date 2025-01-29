

*********************
dysh for GBTIDL Users
*********************


Some key differences
====================

* No longer tied to a single input and output file, can load as many SDFITS files as RAM allows

* No longer tied to a finite number of data containers, you can have as many ScanBlock/Spectrum objects as RAM allows

* Python allows for well-trackable local and global variables





Rough equivalents
=================

In the tables below, it is assumed you have executed the following commands in Python:

`from dysh.fits.gbtfitsload import GBTFITSLoad`




File I/O and Metadata Operations
--------------------------------


 .. csv-table::
    :file: files/FileIO.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10



Calibration and Data Retrieval
------------------------------


 .. csv-table::
    :file: files/Calibration.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10


Spectrum Operations
-------------------


 .. csv-table::
    :file: files/Spectrum_ops.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10
