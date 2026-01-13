.. _scanblocks:

*************************
Data Structures Explained
*************************

GBTFITSLoad 
===========

The main data object in dysh is 
`~dysh.fits.gbtfitsload.GBTFITSLoad`
which represents one or more :ref:`SDFITS files. <sdfits-explanation>`
in the :ref:`GBT SDFITS format <sdfits-reference>`.  (A more generic
SDFITS loader is  `~dysh.fits.gbtfitsload.SDFITSLoad`).  
The data are represented is as one contiguous block even if multiple
files were loaded.  

.. code:: Python

   sdfits = GBTFITSLoad(path)
   sdfits.summary()  # A tabular summary of all the data

where ``path`` is a `path.Path` to a directory containing SDFITS files or a single SDFITS file.  

Once loaded, the columns of the binary table, except the DATA and FLAGS columns, are stored as a `~pandas.DataFrame`, and are accessible via the Python `[]` notation, e.g.:

.. code:: Python

   sdfits["OBJECT"]      # The entire OBJECT column
   set(sdfits["OBJECT"]) # The unique set of OBJECTS
   sdfits.selection         # The entire DataFrame  

This mechanism can be used to :ref:`select data <data-selection>` for calibration.

Although not in the  `~pandas.DataFrame`, The DATA column of a `~dysh.fits.gbtfitsload.GBTFITSLoad` is directly accessible as a `~numpy.ndarray` with ``sdfits["DATA"]``.  More commonly, to look at the raw data you would access a single integration (row) as a numpy array or as a `~dysh.spectrum.Spectrum`, which would include the metadata.

.. code:: Python

   sdfits.rawspectrum(10) #  get data array for row 10 
   sdfits.getspec(10)     #  get a Spectrum for row 10 data and metadata


`~dysh.fits.gbtfitsload.GBTOnline`
`~dysh.fits.gbtfitsload.GBTOffline`

Spectrum
========

Scan
====

ScanBlock
=========

A `~dysh.spectra.scan.ScanBlock` is a `list` of
`~dysh.spectra.scan.ScanBase` objects.  As illustrated
by the diagram below, it serves to hold together a series
of scans.  One of its purposes is to facilitate working with
groups of scans, as it allows time averaging them together (using
:py:meth:`ScanBlock.timeaverage <dysh.spectra.scan.ScanBlock.timeaverage>`),
or subtracting a baseline from the integrations in a group
of scans (using :py:meth:`ScanBlock.subtract_baseline
<dysh.spectra.scan.ScanBlock.subtract_baseline>`).

.. mermaid::
    :caption: The contents of a `~dysh.spectra.scan.ScanBlock` are a series of `~dysh.spectra.scan.ScanBase` objects. `~dysh.spectra.scan.ScanBlock` are the return of the calibration routines (e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad.getps`, `~dysh.fits.gbtfitsload.GBTFITSLoad.getfs` or `~dysh.fits.gbtfitsload.GBTFITSLoad.gettp`)

    flowchart TD

    subgraph newLines[GBTFITSLoad
        50 Scans
        Position Switching
        Dual Linear Polarization
        1 Beam
        4 Frequency Windows
        100 integrations each
        ]




    end


    newLines -- getps( scan=45, plnum=1, ifnum=0, fdnum=0 ) --> ScanBlock1
    newLines -- gettp( scan=[17,18,19], intnum=np.r_[50:100], ifnum=2, plnum=0,fdnum=0 ) --> ScanBlock2



    subgraph ScanBlock1[ScanBlock]
            psscan["spectra.scan.PSScan<br />scans = 44,45<br />plnum = 1<br />ifnum = 0<br />fdnum = 0<br />intnum = (0,100)"]
        end
    subgraph ScanBlock2[ScanBlock]
            tpscan1["spectra.scan.TPScan<br />scan=17<br />plnum = 0<br />ifnum = 2<br />fdnum = 0 <br />intnum=(50,100)"]
            tpscan2["spectra.scan.TPScan<br />scan=18<br />plnum = 0<br />ifnum = 2<br />fdnum = 0 <br />intnum=(50,100)"]
            tpscan3["spectra.scan.TPScan<br />scan=19<br />plnum = 0<br />ifnum = 2<br />fdnum = 0 <br />intnum=(50,100)"]

        end

    ScanBlock1[Scan Block] -- timeaverage() --->spectrum1[Spectrum]
    ScanBlock2[Scan Block] -- timeaverage() --->spectrum2[Spectrum]


The return object of all calibration methods (e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad.getps`, `~dysh.fits.gbtfitsload.GBTFITSLoad.getfs` or `~dysh.fits.gbtfitsload.GBTFITSLoad.gettp`) is a `~dysh.spectra.scan.ScanBlock`.
Since a `~dysh.spectra.scan.ScanBlock` is a `list`, it is possible to append to a `~dysh.spectra.scan.ScanBlock`.
This can be useful when working with different observing procedures or for combining polarizations and/or spectral windows (IFs).
For example, you can do the following:

.. code:: python3

    # Calibrate plnum 1 and 2 data, then average.
    sb0 = sdfits.getps(scan=1, ifnum=0, plnum=0, fdnum=0)
    sb1 = sdfits.getps(scan=1, ifnum=0, plnum=1, fdnum=0)
    sb0.extend(sb1)
    pol_avg_spectrum = sb0.timeaverage()



