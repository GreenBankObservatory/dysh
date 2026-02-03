
.. _usersguide-loading-and-examining-data_:

**************************
Loading and Examining data
**************************

.. _usersguide-loading-sdfits:

The basics classes for reading in SDFITS files were explained 
in :ref:`the previous section <usersguide-gbtfitsload>`.  Here we describe some of 
their keywords and how to inspect data  after loading it.

Header Data Units
-----------------

By default `~dysh.fits.gbtfitsload.GBTFITSLoad` will load all `Header Data Units (HDUs) <https://docs.astropy.org/en/stable/io/fits/api/hdus.html>`_ 
from the input file(s).  In SDFITS format, the HDUs are `binary tables <https://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.BinTableHDU>`_.
For files with multiple HDUs, you can load a specific HDU with the `hdu` keyword, which can speed up data loading if you are only interested in one binary table.

.. code:: Python

    # Load the second binary table from the input FITS files.  
    # Note binary tables start at index 1, since HDU 0 is the 
    # Primary HDU contains only metadata.
    sdfits = GBTFITSLoad("/path/to/data",hdu=2)
    sdfits.info()

You can also limit 
Flags and Flag Files
--------------------

Data may come with flag files (with extension ".flag").  By default dysh does 
not read these files because primarily the contain VEGAS spur channels which are more quickly
flagged algorithmically based on information in the SDFITS header.  You can control 
flagging on input with the `skipflags` and `flag_vegas` keywords, both of which default to True.

.. code:: Python

   # Load a single SDFITS file. Read in .flag file if it exists.
   # VEGAS spurs are still flagged algorithmically.
   sdfits = GBTFITSLoad("/path/to/mydata.fits", skipflags=False)
   
   # Load multiple SDFITS files from a given directory. 
   # Do not read in any .flag files and do not flag VEGAS spurs.
   sdfits = GBTFITSLoad("/path/to/data/", flag_vegas=False)

.. tip::

    See the `flagging chapter <flagging.html>`_ for more details on setting and applying flags.
 
One you have loaded data, you can see what files were read in:

.. code:: Python

   # Print out the files that were loaded
   sdfits.filenames()
   
   # or as Paths instead of strings
   sdfits
   
The most basic description of the data is from `~dysh.fits.gbtfitsload.GBTFITSLoad.summary`, the output of which is customizable.

.. code:: Python

   # The default, compact view
   sdfits.summary()

   # Show every record for a set of scans
   sdfits.summary(scan=[19,20,21],verbose=True)

   # List specific columns
   sdfits.summary(columns=["OBJECT","LST","DATE-OBS"]))

   # List the default columns plus additional ones
   sdfits.summary(add_columns=["TAMBIENT","VELDEF"])

Data in the individual SDFITS files are accessible via the
`~dysh.fits.gbtfitsload.GBTFITSLoad.sdf` attribute, a list
of `~dysh.fits.sdfitsload.SDFITSLoad` with length equal to the number of files.

.. code:: Python

   sdfits.sdf[0]
   # get array of data from first FITS file, row 3
   array = sdfits.sdf[0].rawspectrum(3)

.. _usersguide-examining-metadata:

Examining and Setting Metadata
------------------------------

Metadata are any columns in the SDFITS file other than DATA and FLAGS.
A powerful mechanism for examining and modifying the metadata columns is the `[] operator <https://docs.python.org/3/reference/datamodel.html#object.__getitem__>`_. This is used to access any column of the metadata.

.. code:: Python

   # Return the names of all the columns (except DATA and FLAGS)
   sdfits.columns

   # Examine individual columns, note column names are case-insensitive
   sdfits["obstype"]
   sdfits["tunit7"]
   sdfits["data"].shape

   # Examine multiple columns at once
   mylist = ["object","backend","intnum"]
   sdfits[mylist]

   # The unique set of values in a column
   sdfits.udata["backend"]


Assignment also works. Assigned values will be used in any subsequent calibration commands. The underlying FITS files are not affected unless you actually overwrite them with ``sdfits.write(overwrite=True)``.

.. code:: Python

   # Assign one value to all rows
   sdfits["TCAL"] = 1.5

   # Assign an array with lengths equal to number of rows
   # (a silly example)
   tcal = np.arange(sdfits.stats()['nrows'])
   sdfits["TCAL"] = tcal

New columns can be created by assignment as well, either one value assigned to all rows, or with an array with length equal to the number of rows.

.. code:: Python

   sdfits["PI"] = np.pi
   sdfits["RANDOM"] = np.random.rand(sdfits.stats()['nrows'])

Columns can be renamed:

.. code:: Python

   sdfits.rename_column("PI","PIE")


.. tip::

    For more information on setting metadata, see `Metadata Management <https://dysh.readthedocs.io/en/latest/how-tos/examples/metadata_management.html>`_

.. _usersguide-examining-raw-spectral-data:

Examining the Raw Spectral Data
--------------------------------

GBTFITSLoad, GBTOnline, and GBTOffline have two methods to look at the uncalibrated integration, `~dysh.fits.gbtfitsload.GBTFITSLoad.rawspectrum` returns a `~numpy.ndarray` for a given integration, while `~dysh.fits.gbtfitsload.GBTFITSLoad.getspec` returns the data in a `~dysh.spectra.spectrum.Spectrum` object.
By default these retrieve the the record from the first SDFITS file; other files can be accessed with the ``fitsindex`` parameter (equivalent to, e.g.,  ``sdfits.sdf[2].rawspectrum()``.

.. code:: Python

   array = sdfits.rawspectrum(10) #  get data array for row 10 from the first SDFITS file
   array = sdfits.rawspectrum(10, fitsindex=2) #  get data array for row 10 from the third SDFITS file
   spectrum = sdfits.getspec(10)  #  get a Spectrum for row 10 data and metadata
   spectrum.plot()                #  Spectrum objects can always be plotted.

The full raw data array for any binary table from one of the underlying SDFITS files can be retrieved 
with `~dysh.fits.gbtfitsload.GBTFITSLoad.rawspectra`.

.. code:: Python

    # get the data array from the second binary table of the third SDFITS file
    array = sdfits.rawspectra(bintable=1, fitsindex=2) 
     
.. warning::
   
   `~dysh.fits.gbtfitsload.GBTFITSLoad.rawspectrum` and `~dysh.fits.gbtfitsload.GBTFITSLoad.rawspectra` return references to the actual binary table data.  
   If you alter the result, you alter the data!  It is safer to use the "DATA" keyword.

If there is a single binary table, the entire raw data array can be retrieved using the "DATA" keyword.

   allthedata = sdfits["DATA"]

Note this is copy of the data not a reference to it, modifying ``allthedata`` does not affect the binary table data.


.. _usersguide-adding-comments:

Adding Comments
---------------

Users can add COMMENT or HISTORY cards to the the `~dysh.fits.gbfitsload.GBTFITSLoad` object  with
`~dysh.fits.gbtfitsload.GBTFITSLoad.add_comment`
and
`~dysh.fits.gbtfitsload.GBTFITSLoad.add_history`. These get propagated to ScanBlocks and Spectrum objects during processing and written to the output file(s) during `~dysh.fits.gbtfitsload.GBTFITSLoad.sdfits.write`.
