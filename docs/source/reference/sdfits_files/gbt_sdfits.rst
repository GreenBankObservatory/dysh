.. _sdfits-reference:

****************
GBT SDFITS Files
****************

Convention
==========
The single-dish FITS (SDFITS) convention is used for observer-facing GBT data. The GBT records raw data to FITS files, and then the SDFITS filler populates SDFITS files with information the observer needs. GBT SDFITS files contain 2 or more Header-Data Units (HDUs). Note that keys sometimes appear to be truncated because they can only have up to 8 letters.

HDU 0 (PRIMARY)
---------------

.. _primary-sdfits-header:

Header
^^^^^^

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - DATE
     - UTC Timestamp
     - Date and time the file was created
   * - ORIGIN
     - "NRAO Green Bank"
     - Origin of the observation
   * - TELESCOP
     - "NRAO_GBT"
     - The telescope used for the observation
   * - INSTRUME
     - "VEGAS"
     - The backend used for the observation
   * - GBTMCVER
     - "sdfits ver#.##"
     - Version of SDFITS used to make this file (current: 1.22)
   * - FITSVER
     - "#.#"
     - FITS definition version

HDU 1 (SINGLE DISH)
-------------------

.. _singledish-sdfits-header:

Header
^^^^^^

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - XTENSION
     - "BINTABLE"
     - binary table extension
   * - BITPIX
     - 8
     - 8-bit bytes
   * - NAXIS
     - 2
     - 2-dimensional binary table
   * - NAXIS1
     - integer
     - width of table in bytes
   * - NAXIS2
     - integer
     - number of rows in table
   * - PCOUNT
     - integer
     - size of special data area
   * - GCOUNT
     - 1
     - one data group (required keyword)
   * - TFIELDS
     - integer
     - number of fields in each row
   * - TELESCOP
     - "NRAO_GBT"
     - The telescope used for the observation
   * - PROJID
     - string
     - project identifier
   * - Backend
     - "VEGAS"
     - backend device
   * - SITELONG
     - -7.983983E+01
     - E. longitude of intersection of the az/el axes
   * - SITELAT
     - 3.843312E+01
     - N. latitude of intersection of the az/el axes
   * - SITEELEV
     - 8.245510E+02
     - height of the intersection of az/el axes
   * - EXTNAME
     - "SINGLE DISH"
     - The name of this binary table extension

Data
^^^^

These are the fields (columns) of the BINTABLE.

.. list-table::
   :widths: 20 20 20 40
   :header-rows: 1

   * - Key
     - Unit
     - Type
     - Definition
   * - OBJECT
     -
     - 32A
     - Name of source observed
   * - BANDWID
     - Hz
     - 1D
     - bandwidth
   * - DATE-OBS
     -
     - 22A
     - Start date and time of the integration
   * - DURATION
     - s
     - 1D
     - Total integration duration in seconds
   * - EXPOSURE
     - s
     - 1D
     - Effective int time (excludes blanking) in seconds
   * - TSYS
     - K
     - 1D
     - System temperature in Kelvin
   * - DATA
     -
     - 131072E
     - Actual data (needs to be field 7)
   * - TDIM7
     -
     - 16A
     - Data dimensions of the array in field 7
   * - TUNIT7
     -
     - 6A
     -
   * - CTYPE1
     - Hz
     - 8A
     - First data axis is frequency-like
   * - CRVAL1
     - Hz
     - 1D
     -
   * - CRPIX1
     - index
     - 1D
     -
   * - CDELT1
     - Hz
     - 1D
     -
   * - CTYPE2
     -
     - 4A
     - Second axis is longitude-like
   * - CRVAL2
     - deg
     - 1D
     -
   * - CTYPE3
     -
     - 4A
     - Third axis is latitude-like
   * - CRVAL3
     - deg
     - 1D
     -
   * - CTYPE4
     -
     -
     - Fourth axis is Stokes
   * - CRVAL4
     -
     - 1I
     -
   * - OBSERVER
     -
     - 32A
     - Name of observer(s)
   * - OBSID
     -
     - 32A
     - Observation description
   * - SCAN
     -
     - 1J
     - Scan number
   * - OBSMODE
     -
     - 32A
     - Observing mode
   * - FRONTEND
     -
     - 16A
     - Frontend device
   * - TCAL
     - K
     - 1E
     - Calibration temperature
   * - VELDEF
     -
     - 8A
     - Velocity definition and frame
   * - VFRAME
     - m/s
     - 1D
     - Radial velocity of the reference frame
   * - RVSYS
     - m/s
     - 1D
     - Radial velocity, Vsource - Vtelescope
   * - OBSFREQ
     - Hz
     - 1D
     - Observed center frequency
   * - LST
     - s
     - 1D
     - LST at midpoint of integration/scan
   * - AZIMUTH
     - deg
     - 1D
     - Azimuth
   * - ELEVATIO
     - deg
     - 1D
     - Elevation
   * - TAMBIENT
     - K
     - 1D
     - Ambient temperature
   * - PRESSURE
     - mmHg
     - 1D
     - Ambient pressure
   * - HUMIDITY
     -
     - 1D
     - Relative humidity
   * - RESTFREQ
     - Hz
     - 1D
     - Rest frequency at band center
   * - FREQRES
     - Hz
     - 1D
     - Frequency resolution
   * - EQUINOX
     -
     - 1D
     - Equinox of selected coordinate reference frame
   * - RADESYS
     -
     - 8A
     - Equatorial coordinate system name
   * - TRGTLONG
     - deg
     - 1D
     - Target longitude in coord. ref. frame
   * - TRGTLAT
     - deg
     - 1D
     - Target latitude in coord. ref. frame
   * - SAMPLER
     -
     - 12A
     - Sampler description (e.g., "A1" or "A1xA3")
   * - FEED
     -
     - 1I
     - (signal) feed number
   * - SRFEED
     -
     - 1I
     - Reference feed number
   * - FEEDXOFF
     - deg
     - 1D
     - Feed XEL offset
   * - FEEDEOFF
     - deg
     - 1D
     - Feed EL offset
   * - SUBREF_STATE
     -
     - 1I
     - Subreflector state (1,0,-1) - 0=moving
   * - SIDEBAND
     -
     - 1A
     - Resulting sideband ('U'pper or 'L'ower)
   * - PROCSEQN
     -
     - 1I
     - Scan sequence number
   * - PROCSIZE
     -
     - 1I
     - Number of scans in procedure
   * - PROCSCAN
     -
     - 16A
     - Scan's role in the procedure
   * - PROCTYPE
     -
     - 16A
     - Type of the procedure
   * - LASTON
     -
     - 1J
     - Last 'on' for position switching
   * - LASTOFF
     -
     - 1J
     - Last 'off' for position switching
   * - TIMESTAMP
     - UTC
     - 22A
     - Date and time of scan start
   * - QD_XEL
     - deg
     - 1D
     - QuadrantDetector cross-elevation offset
   * - QD_EL
     - deg
     - 1D
     - QuadrantDetector elevation offset
   * - QD_BAD
     -
     - 1I
     - QuadrantDetector flag: 0=good,1=bad
   * - QD_METHOD
     -
     - 1A
     - Quad. Det. method A,B,C. Blank indicates none.
   * - VELOCITY
     - m/s
     - 1D
     - line velocity in rest frame
   * - ZEROCHAN
     -
     - 1E
     - Zero channel
   * - DOPFREQ
     - Hz
     - 1D
     - Doppler tracked frequency
   * - ADCSAMPF
     -
     - 1D
     - VEGAS ADC sampler frequency
   * - VSPDELT
     -
     - 1D
     - Channel increment between adjacent VEGAS spurs
   * - VSPRVAL
     -
     - 1D
     - VEGAS spur number at VSPRPIX
   * - VSPRPIX
     -
     - 1D
     - Channel number of VEGAS spur VSPRVAL
   * - SIG
     -
     - 1A
     - Signal is true, reference is false
   * - CAL
     -
     - 1A
     - Cal ON is true, cal OFF is false
   * - CALTYPE
     -
     - 8A
     - LOW or HIGH, may eventually be other types
   * - TWARM
     - K or C
     - 1E
     - 4mm RX ambient load temp (K) or Argus vane temperature (C)
   * - TCOLD
     - K
     - 1E
     - 4mm RX cold load temp (K)
   * - CALPOSITION
     -
     - 16A
     - 4mm RX table position
   * - IFNUM
     -
     - 1I
     - Spectral window (IF) number
   * - PLNUM
     -
     - 1I
     - Polarization number
   * - FDNUM
     -
     - 1I
     - Feed number

Index Files
===========

An index file is a tabular ASCII representation of a subset of the SDFITS binary table column (except the DATA column) written by GBTIDL.  It is used as a faster way to some load binary table data without having to read the FITS file(s). Its structure is described below.

.. list-table:: Header
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - created
     - Day Month DD HH:MM:SS YYYY
     - Date file was created
   * - last_modified
     - Day Month DD HH:MM:SS YYYY
     - Date file was last modified
   * - version
     -
     -
   * - observer
     -
     -
   * - backend
     -
     -
   * - tcal_rx_table
     -
     -
   * - sprotect
     -
     -
   * - created_by
     - "index_writer"
     - Method that created this file

.. list-table:: Rows
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - INDEX
     -
     -
   * - PROJECT
     -
     -
   * - FILE EXT
     -
     -
   * - ROW
     -
     -
   * - SOURCE
     -
     -
   * - PROCEDURE
     -
     -
   * - OBSID
     -
     -
   * - E2ESC
     -
     -
   * - PROCS
     -
     -
   * - SCAN
     -
     -
   * - POL
     -
     -
   * - PLNUM
     -
     -
   * - IFNUM
     -
     -
   * - FEED
     -
     -
   * - FDNUM
     -
     -
   * - INT
     -
     -
   * - NUMCHN
     -
     -
   * - SIG
     -
     -
   * - CAL
     -
     -
   * - SAMPLER
     -
     -
   * - AZIMUTH
     -
     -
   * - ELEVATION
     -
     -
   * - LONGITUDE
     -
     -
   * - LATITUDE
     -
     -
   * - TRGTLONG
     -
     -
   * - TRGTLAT
     -
     -
   * - SUB
     -
     -
   * - LST
     -
     -
   * - CENTFREQ
     -
     -
   * - RESTFREQ
     -
     -
   * - VELOCITY
     -
     -
   * - FREQINT
     -
     -
   * - FREQRES
     -
     -
   * - DATEOBS
     -
     -
   * - TIMESTAMP
     -
     -
   * - BANDWIDTH
     -
     -
   * - EXPOSURE
     -
     -
   * - TSYS
     -
     -
   * - NSAVE
     -
     -
   * - PROCSCAN
     -
     -
   * - PROCTYPE
     -
     -
   * - WCALPOS
     -
     -


Flag Files
==========

Flag files indicate the data that should be ignored. For example, these flags can include the locations of VEGAS spurs. `GBTIDL` sometimes auto-masks data that is flagged in these files immediately upon start.  dysh will read flag files in `~dysh.fits.gbtfitsload.GBTFITSLoad` is instantiated with `skip_flags=False`. dysh calculates the location of VEGAS spurs and flags them (regardless of the `skip_flags` value). This is controlled separately by the `flag_vegas` parameter.

.. list-table:: Header
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - created
     - Day Month DD HH:MM:SS YYYY
     - Date file was created
   * - version
     - 1.0
     - Version of ?
   * - created_by
     - sdfits
     - Created by the SDFITS filler

.. list-table:: Flags
   :widths: 25 25 50
   :header-rows: 1

   * - Key
     - Value
     - Definition
   * - RECNUM
     - integer or "*"
     -
   * - SCAN
     - integer or "*"
     - Scan number
   * - INTNUM
     - integer or "*"
     - Integration number
   * - PLNUM
     - integer or "*"
     - Polarization number
   * - IFNUM
     - integer or "*"
     - Spectral window (IF) number
   * - FDNUM
     - integer or "*"
     - Feed number
   * - BCHAN
     - list of integers
     -
   * - ECHAN
     - list of integers
     -
   * - IDSTRING
     - "VEGAS_SPUR"
     - Type of flag


Other Resources
===============
The full SDFITS documentation for GBO can be found on `the GBT SDFITS Project Wiki <https://safe.nrao.edu/wiki/bin/view/GB/Data/Sdfits>`_. However, this page is out of date and requires a login to view.
