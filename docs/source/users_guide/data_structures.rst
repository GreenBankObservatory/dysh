
*************************
Data Structures Explained
*************************

.. _usersguide-gbtfitsload:

==========

GBTFITSLoad
===========

The main data object in dysh is
`~dysh.fits.gbtfitsload.GBTFITSLoad`
which represents one or more :ref:`SDFITS files. <sdfits-explanation>`
in the :ref:`GBT SDFITS format <sdfits-reference>`.  (A more generic
SDFITS loader is  `~dysh.fits.sdfitsload.SDFITSLoad`, described below).
The data are represented as one contiguous block even if multiple
files were loaded.

.. code:: Python

   sdfits = GBTFITSLoad(path)
   sdfits.summary()   # A tabular summary of all the data

where ``path`` is a `path.Path` to a directory containing SDFITS files or a single SDFITS file.

Once loaded, the columns of the binary table, except the DATA and FLAGS columns, are stored as a `~pandas.DataFrame`, and are accessible via the Python `[]` notation, e.g.:

.. code:: Python

   sdfits["OBJECT"]         # The entire OBJECT column
   sdfits.udata("OBJECT")   # The unique set of OBJECTS
   sdfits.selection         # The entire DataFrame

This mechanism can be used to `select data <selection.html>`_ for calibration.

Although not in the  `~pandas.DataFrame`, The DATA column of a `~dysh.fits.gbtfitsload.GBTFITSLoad` is directly accessible as a `~numpy.ndarray` with ``sdfits["DATA"]``.  More commonly, to look at the raw data you would access a single integration (row) as a numpy array or as a `~dysh.spectrum.Spectrum`, which would include the metadata.

.. code:: Python

   array = sdfits.rawspectrum(10) #  get data array for row 10
   spectrum = sdfits.getspec(10)  #  get a Spectrum for row 10 data and metadata


GBTOnline and GBTOffline
========================

For users at GBO, `~dysh.fits.gbtfitsload.GBTOnline` connects directly to the SDFITS file(s) currently being observed/written at the GBT.  It functions exactly like  `~dysh.fits.gbtfitsload.GBTFITSLoad` and updates its contents automatically as new data are written.

.. code:: Python

   sdfits = GBTOnline()
   sdfits.summary()   # A tabular summary of the current data


`~dysh.fits.gbtfitsload.GBTOffline` connects to a given project on disk, by default in ``/home/sdfits`` (where project data live a GBO).  You can set the SDFITS_DATA environment variable to point to a different location.

.. code:: Python

   # load data from project id AGBT_22A_325_33
   sdfits = GBTOffline('AGBT_22A_325_33')
   sdfits.summary()

.. _usersguide-sdfitsload:

SDFITSLoad
==========

A generic class for loading single SDFITS files is `~dysh.fits.sdfitsload.SDFITSLoad`. 
(Underneath the hood, `~dysh.fits.gbtfitsload.GBTFITSLoad` contains one `~dysh.fits.sdfitsload.SDFITSLoad` for each file within).
GBT users will not likely instantiate one of these objects directly. However, they can
be used to read in a SDFITS file that came from a different telescope.  It has basic methods to inspect data, such as
`~dysh.fits.sdfitsload.SDFITSLoad.summary`
`~dysh.fits.sdfitsload.SDFITSLoad.info`, 
`~dysh.fits.sdfitsload.SDFITSLoad.getspec`,
`~dysh.fits.sdfitsload.SDFITSLoad.rawspectrum`,
`~dysh.fits.sdfitsload.SDFITSLoad.rawspectra`,  as well as the `[]` style accessor for column data.

.. code:: Python

    from dysh.fits import SDFITSLoad

    sdfits = SDFITSLoad('AGBT05B_047_01.raw.acs.fits')
    sdfits.summary()
    
    File:    AGBT05B_047_01.raw.acs.fits
    i= 0
    HDU       1
    BINTABLE: 352 rows x 70 cols with 32768 chans
    Selected  352/352 rows
    Sources:  ['NGC5291']
    RESTFREQ: [1.420405] GHz
    Scans:    [51, 52, 53, 54, 55, 56, 57, 58]
    Npol:     2
    Nint:     176
    
    sdf.info()
    
    Filename: AGBT05B_047_01.raw.acs.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU      12   ()      
    1  SINGLE DISH    1 BinTableHDU    229   352R x 70C   ['32A', '1D', '22A', '1D', '1D', '1D', '32768E', '16A', '6A', '8A', '1D', '1D', '1D', '4A', '1D', '4A', '1D', '1I', '32A', '32A', '1J', '32A', '16A', '1E', '8A', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '8A', '1D', '1D', '12A', '1I', '1I', '1D', '1D', '1I', '1A', '1I', '1I', '16A', '16A', '1J', '1J', '22A', '1D', '1D', '1I', '1A', '1D', '1E', '1D', '1A', '1A', '8A', '1E', '1E', '16A', '1I', '1I', '1I']     
 
    sdfits["IFNUM"]

    0      0
    1      0
    2      0
    3      0
    4      0
          ..
    347    0
    348    0
    349    0
    350    0
    351    0
    Name: IFNUM, Length: 352, dtype: int16    

.. warning::

   Be careful accessing the data of an `~dysh.fits.sdfitsload.SDFITSLoad` with the "DATA" key! It is a reference
   to the underlying data, not a copy. You could overwrite the data with `sdfits["DATA"] = ...`
    
.. _usersguide-spectrum:

Spectrum
========

A `~dysh.spectra.spectrum.Spectrum` is a container representing a
spectrum and its attributes, with data in a brightness unit (e.g.,
`astropy.units.ct`, `astropy.units.K`, `astropy.units.Jy`) and a
spectral axis in frequency or velocity units.  The data are accesible as a (unitless) `~numpy.ndarray` (`~specutils.Spectrum.data`) or as a
`~astropy.units.quantity.Quantity` (`~dysh.spectra.spectrum.Spectrum.flux`).  Spectrum objects have a `~dysh.spectra.spectrum.Spectrum.mask` array which can be set with
`flagging operations <https://dysh.readthedocs.io/en/latest/how-tos/examples/flagging.html>`_. It is based on specutils
`~specutils.Spectrum` class.  It supports most common velocity reference
frames supported by `astropy.coordinates.BaseCoordinateFrame` ('itrs','topo', 'lsrk','icrs','hcrs', etc).  'topo' is a synonym for 'itrs'.
Spectra can be displayed with  `~dysh.spectra.spectrum.Spectrum.plot`, which opens up an interactive plot window.
The frequency/velocity locations of spectral lines within the spectral window can be listed using `~dysh.spectra.spectrum.Spectrum.query_lines`.

Standard operations such as `~dysh.spectra.spectrum.Spectrum.baseline` removal,  `~dysh.spectra.spectrum.Spectrum.smooth`, and `~dysh.spectra.spectrum.Spectrum.average` are supported, as well as analysis functions like  `~dysh.spectra.spectrum.Spectrum.stats`,  `~dysh.spectra.spectrum.Spectrum.roll` ,  `~dysh.spectra.spectrum.Spectrum.radiometer`, `~dysh.spectra.spectrum.Spectrum.normalness`, and  `~dysh.spectra.spectrum.Spectrum.cog` (:ref:`Curve of Growth <cog>`).  Spectrum arithmetic is supported with common operators, e.g.:

.. code:: Python

    s1 = Spectrum.fake_spectrum(nchan=2048)   # default unit is K
    s2 = Spectrum.fake_spectrum(nchan=2048)
    s3 = s1+s2      # Sum of two spectra
    ratio = s2/s1   # Unitless ratio

`~dysh.spectra.spectrum.Spectrum` objects are returned by operations like `~dysh.spectra.scan.ScanBase.timeaverage` and `~dysh.fits.gbtfitsload.GBTFITSLoad.getspec`.
A useful method for creating dummy spectra is `~dysh.spectra.spectrum.Spectrum.fake_spectrum`.


.. _usersguide-scan:

Scan
====

Scans are represented by subclasses of `~dysh.spectra.scan.ScanBase` specific to their observation type.  For instance
`~dysh.spectra.scan.PSScan` is position-switched,
`~dysh.spectra.scan.FSScan` is frequency-switched,
`~dysh.spectra.scan.TPScan` is total power,
`~dysh.spectra.scan.NodScan` is nodding.
These are created by the corresponding calibration routine, e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad.getps`,  `~dysh.fits.gbtfitsload.GBTFITSLoad.getfs`, `~dysh.fits.gbtfitsload.GBTFITSLoad.gettp`
and collected into a common :ref:`usersguide-scanblock`.  Users generally will not interact with individual Scan objects, but rather collectively through the ScanBlock interface.  However, certain attributes and functions are common to all `Scan classes <https://dysh.readthedocs.io/en/latest/reference/modules/dysh.spectra.html#dysh.spectra.scan.ScanBase>`_ are helpful for understanding the data.

.. code:: Python

   # Load data from project id AGBT_22A_325_33
   sdfits = GBTFITSLoad(path)

   # Get a ScanBlock containing some PSScans calibrated to Janskys
   sb = sdfits.getps(scan=[23,25,27], ifnum=0,plnum=0,fdnum=0, units='flux', zenith_opacity=0.05)

   # Look at the scale factors and weights of all the Scans in the ScanBlock
   for s in sb:
       print(f"Scale factors for PSScan {s.scan} = {s.tscale_fac}")
       print(f"Weights for PSScan {s.scan} = {s.weights}")

   # Get the weighted average of the integrations for the first PSScan
   ta = sb[0].timeaverage()
   ta.plot()

   # Examine the integrations in the first PSScan
   for i in range(len(sb[0]):
        sb[0].getspec(i).plot()

   sb[0].add_comment("Integration 12 looks wonky.")

.. _scanblocks:
.. _usersguide-scanblock:

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
