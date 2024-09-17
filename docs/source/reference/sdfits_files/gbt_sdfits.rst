****************
GBT SDFITS Files
****************

Convention
==========
The single-dish FITS (SDFITS) convention is used for observer-facing GBT data. The GBT records raw data to FITS files, and then the SDFITS filler populates SDFITS files with information the observer needs. GBT SDFITS files contain 2 or more Header-Data Units (HDUs). Note that keys sometimes appear to be truncated because they can only have up to 8 letters.

HDU 0 (PRIMARY)
---------------

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
     - name of source observed
   * - BANDWID
     - Hz
     - 1D
     - bandwidth
   * - DATE-OBS
     -
     - 22A
     - date and time of observation start
   * - DURATION
     - s
     - 1D
     - total integration duration in seconds
   * - EXPOSURE
     - s
     - 1D
     - effective int time (excludes blanking) in seconds
   * - TSYS
     - K
     - 1D
     - system temperature in Kelvin
   * - DATA
     -
     - 131072E
     - Actual data
   * - TDIM7
     -
     - 16A
     - data dimensions of the array
   * - TUNIT7
     -
     - 6A
     -
   * - CTYPE1
     - Hz
     - 8A
     - first data axis is frequency-like
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
     - second axis is longitude-like
   * - CRVAL2
     - deg
     - 1D
     -
   * - CTYPE3
     -
     - 4A
     - third axis is latitude-like
   * - CRVAL3
     - deg
     - 1D
     -
   * - CTYPE4
     -
     -
     - fourth axis is Stokes
   * - CRVAL4
     -
     - 1I
     -
   * - OBSERVER
     -
     - 32A
     - name of observer(s)
   * - OBSID
     -
     - 32A
     - observation description
   * - SCAN
     -
     - 1J
     - scan number
   * - OBSMODE
     -
     - 32A
     - observing mode
   * - FRONTEND
     -
     - 16A
     - frontend device
   * - TCAL
     - K
     - 1E
     - calibration temperature
   * - VELDEF
     -
     - 8A
     - velocity definition and frame
   * - VFRAME
     - m/s
     - 1D
     - radial velocity of the reference frame
   * - RVSYS
     - m/s
     - 1D
     - radial velocity, Vsource - Vtelescope
   * - OBSFREQ
     - Hz
     - 1D
     - observed center frequency
   * - LST
     - s
     - 1D
     - LST at midpoint of integration/scan
   * - AZIMUTH
     - deg
     - 1D
     - azimuth
   * - ELEVATIO
     - deg
     - 1D
     - elevation
   * - TAMBIENT
     - K
     - 1D
     - ambient temperature
   * - PRESSURE
     - mmHg
     - 1D
     - ambient pressure
   * - HUMIDITY
     -
     - 1D
     - relative humidity
   * - RESTFREQ
     - Hz
     - 1D
     - rest frequency at band center
   * - FREQRES
     - Hz
     - 1D
     - frequency resolution
   * - EQUINOX
     -
     - 1D
     - equinox of selected coordinate reference frame
   * - RADESYS
     -
     - 8A
     - Equatorial coordinate system name
   * - TRGTLONG
     - deg
     - 1D
     - target longitude in coord. ref. frame
   * - TRGTLAT
     - deg
     - 1D
     - target latitude in coord. ref. frame
   * - SAMPLER
     -
     - 12A
     - sampler description (e.g., "A1" or "A1xA3")
   * - FEED
     -
     - 1I
     - (signal) feed number
   * - SRFEED
     -
     - 1I
     - reference feed number
   * - FEEDXOFF
     - deg
     - 1D
     - feed XEL offset
   * - FEEDEOFF
     - deg
     - 1D
     - feed EL offset
   * - SUBREF_STATE
     -
     - 1I
     - subreflector state (1,0,-1) - 0=moving
   * - SIDEBAND
     -
     - 1A
     - resulting sideband ('U'pper or 'L'ower)
   * - PROCSEQN
     -
     - 1I
     - scan sequence number
   * - PROCSIZE
     -
     - 1I
     - number of scans in procedure
   * - PROCSCAN
     -
     - 16A
     - scan's role in the procedure
   * - PROCTYPE
     -
     - 16A
     - type of the procedure
   * - LASTON
     -
     - 1J
     - last 'on' for position switching
   * - LASTOFF
     -
     - 1J
     - last 'off' for position switching
   * - TIMESTAMP
     - UTC
     - 22A
     - date and time of scan start
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
     - zero channel
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
     - channel increment between adjacent VEGAS spurs
   * - VSPRVAL
     -
     - 1D
     - VEGAS spur number at VSPRPIX
   * - VSPRPIX
     -
     - 1D
     - channel number of VEGAS spur VSPRVAL
   * - SIG
     -
     - 1A
     - signal is true, reference is false
   * - CAL
     -
     - 1A
     - cal ON is true, cal OFF is false
   * - CALTYPE
     -
     - 8A
     - LOW or HIGH, may eventually be other types
   * - TWARM
     - K
     - 1E
     - 4mm RX ambient load temp (K)
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

Flag files indicate the data that should be ignored. For example, these flags can include the locations of VEGAS spurs. `GBTIDL` sometimes auto-masks data that is flagged in these files immediately upon start.

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
