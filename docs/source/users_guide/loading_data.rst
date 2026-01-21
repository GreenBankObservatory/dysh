
.. _usersguide-loading-and-examining-data_:

**************************
Loading and Examining data
**************************

.. _usersguide-loading-sdfits:

The basics of reading in SDFITS files was explained in :ref:`the previous section <usersguide-gbtfitsload>`.  Below are a few more examples.

.. code:: Python

   from dysh.fits import GBTFITSLoad
   import numpy as np

   # Load a single SDFITS file. Don't read the .flag file if it exists.
   sdfits = GBTFITSLoad("/path/to/mydata.fits",skipflags=True)

   # Load all files with .fits extension in a specfic directory
   sdfits = GBTFITSLoad("/path/to/datafiles")
   # Print out the files that were loaded
   sdfits.filenames()
   # or as Paths instead of strings
   sdfits

   # Load all data from a specific project. This assumes
   # /home/sdfits exists or environment variable $SDFITS_DATA is
   # a directory containing projects.
   sdfits = GBTOffline("myproject_id")


The most basic description of the data is from `~dysh.fits.gbtfitsload.GBTFITSLoad.summary`, the output of which is customizable.

.. code:: Python

   # The default, compact view
   sdfits.summary()

   # Show every record
   sdfits.summary(verbose=True)

   # List specific columns
   sdfits.summary(columns=["OBJECT","LST","DATE-OBS"]))

   # List the default columns plus additional ones
   sdfits.summary(add_columns=["TAMBIENT","VELDEF"])

Data in the individual FITS files are accessible via the
`~dysh.fits.gbtfitsload.GBTFITSLoad.sdf` attribute, a list
`~dysh.fits.sdfitsload.SDFITSLoad` with length equal to the number of files.

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


Assignment also works. Assigned values will be used in any subsequent calibration commands. The underlying FITS files are not affected unless you actually overwrite them with `sdfits.write(overwrite=True)`.

.. code:: Python

   # Assign one value to all rows.
   sdfits["TCAL"] = 1.5

   # Assign an array with lengths equal to number of rows.
   # A silly example.
   tcal = np.arange(sdfits.stats()['nrows'])
   sdfits["TCAL"] = tcal

New columns can be created by assignment as well, either one value assigned to all rows, or with an array with length equal to the number of rows.

.. code:: Python

   sdfits["PI"] = np.pi
   sdfits["RANDOM"] = np.random.rand(sdfits.stats()['nrows'])

Columns can be renamed:

.. code:: Python

   sdfits.rename_column("PI","PIE")

For more information on setting metadata, see `Metadata Management <https://dysh.readthedocs.io/en/latest/how-tos/examples/metadata_management.html>`_

.. _usersguide-examining-raw-spectral-data:

Examining the Raw Spectral Data
--------------------------------

GBTFITSLoad, GBTOnline, and GBTOffline have two methods to look at the uncalibrated integration, `~dysh.fits.gbtfitsload.GBTFITSLoad.rawspectrum` returns a `~numpy.ndarray` for a given integration, while `~dysh.fits.gbtfitsload.GBTFITSLoad.getspec` returns the data in a `~dysh.spectra.spectrum.Spectrum` object.
By default these retrieve the the record from the first FITS file; other files can be accessed with the `fitsindex` parameter (equivalent to, e.g.,  `sdfits.sdf[2].rawspectrum()`.

.. code:: Python

   array = sdfits.rawspectrum(10) #  get data array for row 10 from the first FITS file
   array = sdfits.rawspectrum(10, fitsindex=2) #  get data array for row 10 from the third FITS file
   spectrum = sdfits.getspec(10)  #  get a Spectrum for row 10 data and metadata
   spectrum.plot()                #  Spectrum objects can always be plotted.

The entire raw data array can be retrieved using the "DATA" keyword.

   allthedata = sdfits["DATA"]

Note this is copy of the data not a reference to it, modifying `allthedata` does not affect the binary table data.

.. _usersguide-adding-comments:

Adding Comments
---------------

Users can add COMMENT or HISTORY cards to the the `~dysh.fits.gbfitsload.GBTFITSLoad` object  with
`~dysh.fits.gbtfitsload.GBTFITSLoad.add_comment`
and
`~dysh.fits.gbtfitsload.GBTFITSLoad.add_history`. These get propagated to ScanBlocks and Spectrum objects during processing and written to the output file(s) during `~dysh.fits.gbtfitsload.GBTFITSLoad.sdfits.write`.
