.. _glossary:

A dysh glossary
---------------

See also the `GBO Glossary <https://gbtdocs.readthedocs.io/en/latest/glossary.html>`_


.. glossary::


    Argus
      A 16-:term:`pixel` W-band (74-116 GHz) focal plane array in use at the GBT. Named after a mythical figure
      with 100 eyes. See also https://www.gb.nrao.edu/argus/

    band
      A coherent section of channels in frequency space, all with
      the same channel width. Sometimes called an IF band. In radio
      astronomy bands are often referred to by an alphabetical designation,
      e.g. L-band covers 1-2 GHz. A summary
      of the bands commonly used at GBO can be found on
      https://gbtdocs.readthedocs.io/en/latest/references/receivers.html#gregorian-receivers

      See also :term:`ifnum` and :term:`IF`

    bank
      Overloaded term for a **band**, possibly referring to hardware.
      For the GBT a bank is one of VEGAS roaches (or more than one?). The closest
      thing is an IF or VEGAS subband.

    baseline
      Baseline is a generic term usually taken to mean the
      instrumental plus continuum bandpass shape in an observed
      spectrum, or changes in the background level in a continuum
      observation.

    beam
      The footprint of one receiver horn on the sky. Argus has a
      4x4 multi-beam receiver, numbered 0 through 15.

      See also :term:`fdnum`

      See also :term:`multi-beam`

    Beam Switching
      This is a variation on position switching using a receiver
      with multiple beams. The :term:`SIG` and :term:`REF`
      positions on the sky are
      calculated so that the receiver is always pointing at the source. This is most
      useful for point sources. At GBT only implemented for Ka-band.
      See also :term:`Position Switching`

    BINTABLE
      Binary table. In dysh data, BINTABLE is an index running from 0 to N-1,
      where N is the number of binary tables in the SDFITS file.
      see also :term:`FITS`

    blanking
      blanking is a term used in the (VEGAS) correllator, where bad data has been replaced
      with a Not-a-Number value. Not to be confused with the concepts :term:`flagging`
      and :term:`masking` in dysh

    CAL
      see also :term:`SIG` or :term:`REF`
      ON/OFF   SIG/REF 

    caloff
      Signal with no calibration diode in the signal path.

    calon
      Signal with a calibration diode in the signal path.

    Chebyshev
      a type of orthogonal polynomial that is commonly used in
      numerical methods due to its optimal convergence properties and
      connection to the Fourier transform.

    cog - Curve of Growth
      ...
      see also notebook

    DYSH_DATA
      (optional) environment variable pointing to a directory with local copies
      for developers.

      See also :term:`SDFITS_DATA`.

    ECSV
      (Enhanced Character Separated Values) a self-describing ASCII table format popularized by astropy.
      See also https://github.com/astropy/astropy-APEs/blob/main/APE6.rst

    fdnum
      Feed Number. 0, 1, ...
      Also used as the **fdnum=** keyword in getXX()

      See also :term:`beam`

    FITS
      (Flexible Image Transport System): the export format
      for data-cube, although there is also a waterfall cube
      (time-freq-pixel) cube available.

    flagging
      flagging is a non-destructive operation, where data in the
      time-frequency domain is flagged to be skipped.

      Flagging specific to the VEGAS backend, which has bad channels
      also known as 'spurs' at regular channel intervals. VEGAS
      flagging is done automatically by
      :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`.

      See also :term:`masking`


    flag files
      SDFITS files created by GBTIDL can have a separate ASCII flag
      file. By default, :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`
      reads this file and applies the flags therein.

    FWHM
      (Full Width at Half Max): the effective resolution of the
      beam is normally given in :term:`FITS` keywords
      BMAJ, BMIN, and BPA.

    Frequency Switching
      This is a variation on position switching using a receiver
      where the IF is changed.
      See also :term:`Position Switching`

    GBTIDL
      Green Bank Telescope Interactive Data Language. The GBT data
      reduction package written in :term:`IDL` for analyzing GBT spectral line
      data.

    getXX()
      Generic name for the dysh calibration routines, e.g. getps, getfs, getnod etc.

    horn
      Another term used for :term:`beam` or :term:`pixel`.

    IDL
      The Interactive Data Language program, currently of ITT Visual Information Solutions
      but with a long history of owners.

    IF
      Intermediate Frequency, is a frequency to which a carrier wave is shifted as
      an intermediate step in transmission or reception. The terms
      See also :term:`band` and :term:`window` are often used as well, where they
      mean an IF band.

    ifnum
      IF number (0,1,...)
      Also used as the **ifnum=** keyword in getXX().

      See also :term:`band` and :term:`window`

    intnum
      Integration number. 0 being the first.
      Also used as the **intnum=** keyword in getXX()

    KFPA
      K-band Focal Plane Array, a hexagonal set of beams, with a central beam. Covers 18-26 GHz.
      See the `KFPA receiver page <https://gbtdocs.readthedocs.io/en/latest/references/receivers/kfpa.html>`_
      for more details.
   
    masking
      Masking removes or hides the value in the spectrum. As in numpy,
      as mask value of True means the underlying value is not used. In
      dysh masks are set on individual integrations during calibration
      [getXX()]; resultant spectra will have the final mask set in
      Spectrum.mask. See also :term:`flagging`


    metadata
      describes data. Examples for a spectrum are the RA and DEC associated with the spectrum.
      Typically GBT spectra have 70 items in the metadata, implemented as columns in the
      :term:`BINTABLE`
      and accessed via keyword in :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`, e.g., sdf["object"].

      dysh spectra have metadata in Spectrum.meta and Scans in Scan.meta.

    multi-beam
      If an instrument has multiple :term:`beam`s that typically point are different areas in the sky
      (e.g. :term:`Argus` in a 4x4 configuration, and :term:`KFPA` in a 7 beam hexagonal shape).

    Nod or Nodding
      An observing mode where two beams alternatingly look at source and (different) sky.

    Noise Diode
      A device with known effective temperature that is coupled to the
      telescope system to give a measure of system temperature
      (Tsys). When the telescope is pointed on blank sky, the noise
      diode is turned on and then off to determine the off-source
      system temperature. This device is also refered to as the "Cal".
      See also :term:`calon` and :term:`caloff` and 

    OTF Mapping
      On-the-fly mapping: in this procedure the telescope is scanned across the sky to
      sample the emission. The samples are then "gridded" into a map (which is not part
      of dysh). See for example `gbtgridder <https://github.com/GreenBankObservatory/gbtgridder>`_

    pixel
      An overloaded term. Sometimes referred to as the :term:`beam`, but usually interpreted
      in image processing as
      the size of a single (usually square) element in a gridded map (e.g. from an OTF), which
      is commonly referred to as a *picture element*.

    plnum
      Polarization number (0,1,...). Usually 0 and 1, but of course up to 4 values could be present
      for a full Stokes. Averaging the two polarizations will reduce the noise by :math:`sqrt{2}`

      Also used as the **plnum=** keyword in getXX()


    Position Switching
      This is a standard way to obtain spectra by switching
      between a :term:`SIG` and :term:`REF` position on the sky,
      usually using a single beam. For our
      multi-beam receivers see also :term:`Beam Switching`


    Project ID
      A code designating the year and proposal number, e.g. GBT21B-024.  Data associated with
      a project are found in /home/sdfits (or $SDFITS_DATA), with a slight twist of the name.
      In the example this becomes AGBT21B_024.
      See below :ref:`data_org`

    REF
      Reference point, meant to have no signal. See also :term:`CAL`

    Region
      Region or regions of spectrum, use for flagging/masking,baseline subtraction.

    Scan
       A unit of observing, usually in some common mode, with one or more integrations.
       GBT differentiates between different types of scans. Scans are integers,
       starting with 1.


    ScanBlock
      A container for a series of **scan**'s.

      See also :term:`scanblocks`

    SDFITS
      Single Dish **FITS** format, normally used to store
      raw or even calibrated spectra in a FITS binary table (BINTABLE) format.  Each
      row in a BINTABLE has an attached RA,DEC (and other meta-data),
      plus the whole spectrum. This standard was drafted in 1995 (Liszt),
      and has been implemented by many telescopes (Arecibo, FAST, GBT, Parkes, ....),
      albeit with slightly different conventions.  Also to note is that an SDFITS file
      can have more than one BINTABLE extension.

      See also :ref:`sdfits-reference`

    SDFITS_DATA
      (optional) environment variable pointing to a directory where SDFITS
      project directories and files are stored.

    SESSION
      see :ref:`data_org`

    SFL
      Sanson-Flamsteed projection, sometimes used in gridding OTF maps.
      (the GLS - GLobal Sinusoidal is similar to SFL).

    SIG
      signal - see also CAL.

    Spectral Window
      This is closest to what we call a **bank**,
      or **band**, a set of linearly spaced channels.

      See also :term:`ifnum`

    Spectrum
      A coherent section in frequency space, with its own unique meta-data (such as polarization,
      ra, dec, time). Normally the smallest portion of data we can assign. A spectrum is
      defined by its own seting of *(crval, crpix, cdelt)* in a FITS WCS sense.

    SubBeamNod
      Subreflector Beam Nodding. The getXX() is now called `subbeamnod`

    tcal
      Derive the noise diode temperature from observations

    VEGAS
      Versatile GBT Astronomical Spectrometer - https://www.gb.nrao.edu/vegas/

    

..    The velocity of a source using the relativistic definition of the velocity-frequency relationship.

..    The velocity of a source using the optical definition of the velocity-frequency relationship.

..    The velocity of a source using the radio definition of the velocity-frequency relationship.


    waterfall plot
      A plot (or two-dimensional image) that shows time vs. frequency.

    Window
      See :term:`Spectral Window`


.. _data_org:

Data : Project ID / Session
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally projects are assigned a project id, e.g. *AGBT21B_024*, which is
then observed in a number of sessions, numbered starting with 1. The SDFITS data associated
with these are stored under **$SDFITS_DATA**, e.g. for session 5 in this example, this would be
in **$SDFITS_DATA/AGBT21B_024_05/**.

One possible confusion: a project named "GBT21B-024", is labeled "AGBT21B_024" as the
filename prefix for file storage, which is the name that users need for dysh.


.. bands listed alpabetically in the GBO glossary
.. C   4-8 GHz
.. K   18-26
.. Ka  26-40
.. Ku  12-18
.. L   1-2
.. P   300-1000 MHz
.. Q   40-50
.. S   2-4
.. W   75-111
.. X   8-12
