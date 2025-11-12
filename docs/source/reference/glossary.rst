.. _glossary:

A dysh glossary
---------------

See also the `GBO Glossary <https://gbtdocs.readthedocs.io/en/latest/glossary.html>`_

In this glossary we also note overloaded terms.


.. glossary::
    :sorted:

    Argus
      A 16-:term:`pixel` W-band (74-116 GHz) focal plane array in use at the GBT.
      Named after a mythical figure with 100 eyes. See also https://www.gb.nrao.edu/argus/

    band
      A contiguous section of the radio spectrum.
      In radio astronomy bands are often referred to by an alphabetical designation,
      e.g. L-band covers 1-2 GHz.
      A summary of the bands commonly used at GBO can be found on
      https://gbtdocs.readthedocs.io/en/latest/references/receivers.html

    bank
      A VEGAS bank is a single ROACH-2 board with 2 inputs,
      typically for two polarizations of one :term:`IF`.
      A VEGAS bank can produce between one and eight :term:`spectral windows <spectral window>`, which are identified by an :term:`ifnum`.

    baseline
      Baseline is a generic term usually taken to mean the
      instrumental plus continuum bandpass shape in an observed
      spectrum, or changes in the background level in a continuum
      observation.

    beam
      The footprint of one receiver horn on the sky.

      See also :term:`fdnum`, :term:`multi-beam`, and :term:`horn`.

    beam switching
      This is a variation on :term:`position switching` using a receiver
      with multiple beams. The signal and reference positions are chosen so
      that one of the feeds is always pointing at the source. Using this method
      requires that the source be smaller than the on-sky separation between feeds.
      :term:`Nodding <nod>` and :term:`subreflector beam nodding <subbeamnod>` are forms of beam switching.

      See also :term:`position switching`.

    BINTABLE
      Binary table. In dysh, ``bintable`` is an index running from 0 to N-1,
      where N is the number of binary tables in the :term:`SDFITS` file.

      See also :term:`FITS`.

    blanking
      Blanking is the process of replacing data values with a blank.
      For :term:`SDFITS` files, the blanking value is Not-a-Number (NaN).
      For example :term:`VEGAS` will blank data while it is switching states (`Kepley et al. 2015 <https://library.nrao.edu/public/memos/gbt/GBT_288.pdf>`_).
      Not to be confused with the concepts :term:`flagging` and :term:`masking` in dysh.

    CAL
      Column in :term:`SDFITS` files indicating if the noise diode was being fired at a particular time.

    caloff
      Signal with no calibration diode in the signal path (`CAL="F"`).

    calon
      Signal with a calibration diode in the signal path (`CAL="T"`).

    Chebyshev
      A type of orthogonal polynomial that is commonly used in
      numerical methods due to its optimal convergence properties and
      connection to the Fourier transform. One of the options in
      baseline fitting in dysh.

    CoG
      Curve of Growth: integrating the flux from the center of a line outwards.

      See also :ref:`cog` for the dysh implementation.


    DYSH_DATA
      (optional) environment variable pointing to a directory with local copies
      of SDFITS data for developers.

      See also :term:`SDFITS_DATA`.

    ECSV
      Enhanced Character Separated Values: a self-describing ASCII table format popularized by astropy.
      See also https://github.com/astropy/astropy-APEs/blob/main/APE6.rst

    fdnum
      Feed number. An integer, starting at 0, used to identify different :term:`beams <beam>`.
      Also used as the ``fdnum`` keyword in the :term:`getXX()` routines.
      For example, to select the first :term:`beam` one would use ``fdnum=0``.

      See also :term:`beam`.

    FITS
      Flexible Image Transport System: the export format
      for data-cubes and images.

    flagging
      Flagging is a non-destructive operation, where data in the
      time-frequency domain is flagged to be skipped, :term:`masked <masking>` or :term:`blanked <blanking>`.

      Flagging specific to the :term:`VEGAS` backend, which has bad channels
      also known as 'spurs' at regular channel intervals. :term:`VEGAS`
      flagging is done automatically by
      :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`.

      The data are flagged by :class:`~dysh.fits.gbtfitsload.GBTFITSLoad` (or the user).
      Masking is the application of flags using ``apply_flags``.

      See also :term:`masking` and :term:`blanking`.

    flag files
      SDFITS files created by GBTIDL can have a separate ASCII flag
      file. By default, :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`
      reads this file and applies the flags therein.

    FWHM
      Full Width at Half Max.
      A measure of the width of a curve. It reports the width of the
      curve at its half power point. It is commonly used to describe
      the angular resolution of a telescope (also referred to as half
      power beam width, HPBW, in this case), or the width of a
      spectral line.

      The :term:`FITS` keywords BMAJ, BMIN, and BPA  are used for the
      major axis, minor axis, and position angle respectively when referring
      to a spatial beam.


    frequency switching
      This is an observing technique in which a local oscillator is alternated to switch the :term:`IF` into signal and reference states.

    GBTIDL
      Green Bank Telescope Interactive Data Language. The GBT data
      reduction package written in :term:`IDL` for analyzing GBT spectral line
      data.

    getXX()
      Generic name for any of the dysh calibration routines, e.g. getps, getfs, getnod etc.

    horn
      Another term used for :term:`beam` or :term:`pixel`.

    IDL
      The Interactive Data Language program, currently of ITT Visual Information Solutions
      but with a long history of owners.

    IF
      Intermediate Frequency, is a frequency to which a carrier wave is shifted as
      an intermediate step in transmission or reception.
      The term :term:`window` is often used as well, where it means the contiguous range of frequencies being recorded by the spectrometer (e.g., :term:`VEGAS`).

      See also :term:`ifnum`.

    ifnum
      IF number (0,1,...). These are used to identify :term:`spectral windows <spectral window>`.
      Also used as the ``ifnum`` keyword in :term:`getXX()`. For example, to select the first spectral
      window one would use ``ifnum=0``.

      See also :term:`band` and :term:`window`.

    intnum
      Integration number, starting at 0. Also used as the ``intnum`` keyword in the
      :term:`getXX()` routines. For example, to select the first integration one would use ``intnum=0``.

    KFPA
      K-band Focal Plane Array, a hexagonal set of beams, with a central beam. Covers 18-26 GHz.
      See the `KFPA receiver page <https://gbtdocs.readthedocs.io/en/latest/references/receivers/kfpa.html>`_
      for more details.

    masking
      Masking hides the value in the spectrum. As in numpy,
      as mask value of True means the underlying value is not used. In
      dysh masks are set on individual integrations during calibration
      [getXX()]; resultant spectra will have the final mask set in
      ``Spectrum.mask``.

      See also :term:`flagging`.


    metadata
      Describes data. Examples for a spectrum are the right ascencion (RA) and declination (DEC) associated with the spectrum.
      Typically GBT spectra have 70 items in the metadata, implemented as columns in the
      :term:`BINTABLE`
      and accessed via keyword in :class:`~dysh.fits.gbtfitsload.GBTFITSLoad`, e.g., sdf["object"].

      dysh spectra have metadata in `Spectrum.meta` and Scans in `Scan.meta`.

    multi-beam
      If an instrument has multiple :term:`beams <beam>` that typically point to different sky locations
      (e.g. :term:`Argus` in a 4x4 configuration, and :term:`KFPA` in a 7 beam hexagonal shape).

    nod
      An observing mode where two beams alternatingly look at source and (different) sky.

      See also :term:`position switching` and :term:`beam switching`.

    noise diode
      A device with known effective temperature (see :term:`tcal`) that is coupled to the
      receiver to give a measure of system temperature (TSYS).
      When the telescope is pointed on blank sky, the noise diode can be set toalternate between
      on and off states to determine the system temperature.
      This device is also refered to as the "Cal".

      See also :term:`calon` and :term:`caloff`.

    ON/OFF
      In the context of position switching, the on-source (target or signal) and off-source (reference) positions.
      In the context of the noise diode, it being on (firing) or off (not firing).
      The ON/OFF references are an overloaded term for when we refer to the
      :term:`SIG` and :term:`REF` respectively.

    OTF mapping
      On-the-fly mapping: in this procedure the telescope is scanned across the sky to
      sample the emission. The samples are "gridded" on to a map using the tool
      `gbtgridder <https://github.com/GreenBankObservatory/gbtgridder>`_. The gridding
      is not implemented in dysh.

    pixel
      An overloaded term. Sometimes referred to as the :term:`beam`, but usually interpreted
      in image processing as
      the size of a single (usually square) element in a gridded map (e.g. from an OTF), which
      is commonly referred to as a *picture element*.

    plnum
      Polarization number (0,1,...). An integer used to identify the polarization. Usually 0 and 1, but of course up to 4 values could be present
      for full Stokes observations. Averaging two orthogonal polarizations should reduce the noise by :math:`\sqrt{2}`

      Also used as the ``plnum`` keyword in :term:`getXX()`. For example, one would use ``plnum=0`` to select the first polarization.

    position switching
      This is a standard way to obtain spectra by switching
      between a :term:`SIG` and :term:`REF` position on the sky,
      usually using a single beam.

      See also :term:`beam switching` and :ref:`sdmath`

      .. See also :numref:`Equation %s  <eq_sdmath3>`

      .. See also :numref:`eq_sdmath3`

    project code
      A code designating the year and proposal number, e.g. GBT21B-024.  Data associated with
      a project are found in /home/sdfits (or $SDFITS_DATA), with a slight twist of the name.
      In the example this becomes AGBT21B_024.

      See also :ref:`data_org`

    REF
      Reference point, meant to have no signal.

      See also :ref:`sdmath`.

    region
      Region or regions of spectrum, used for flagging/masking, or baseline subtraction.

    scan
       A unit of observing, usually in some common mode, with one or more integrations.
       Scans are referred to as 1-based integers. The observing procedures (e.g., OnOff or Track)
       are commonly referred to as scans (e.g., an OnOff scan or a Track scan).
       A scan can contain multiple integrations.

    ScanBlock
      A container for a series of :term:`scans <scan>`.

      See also :ref:`scanblocks` for more details.

    SDFITS
      Single Dish :term:`FITS` format, normally used to store
      raw or even calibrated spectra in a :term:`FITS` binary table (BINTABLE) format.  Each
      row in a BINTABLE has an attached RA,DEC (and other meta-data),
      plus the whole spectrum. This standard was drafted in 1995 (Liszt),
      and has been implemented by many telescopes (Arecibo, FAST, GBT, Parkes, ....),
      albeit with slightly different conventions.
      An SDFITS file can have more than one BINTABLE extension.

      See also :ref:`sdfits-reference`.

    SDFITS_DATA
      (optional) environment variable pointing to a directory where :term:`SDFITS`
      project directories and files are stored.

    session
      Or Session ID. This is the number (starting at index 01) denoting the observing sessions
      within a given project. It is concatenated with the :term:`project code` to define a unique
      observing session, for example AGBT21B_024_01.

      See also :ref:`data_org`.

    SIG
      Signal.
      In the context of position switching the ON or target position in an ON/OFF observation.
      In the context of frequency switching, the signal state has SIG="T" and the reference state SIG="F".

      See also :ref:`sdmath`.

    spectral window
      Subdivision of a :term:`band` or :term:`bank` to a set of linearly spaced channels in frequency space.
      For GBT data a spectral window is also refered to as an IF.

      See also :term:`ifnum`.

    spectrum
      A coherent section in frequency space, with its own unique meta-data (such as polarization,
      right ascencion, declination, and time). Normally the smallest portion of data we can assign. A spectrum is
      defined by its own seting of *(crval, crpix, cdelt)* in a FITS WCS sense.

    subbeamnod
      Subreflector beam nodding.
      In this observing method the subreflector is used to alternate the signal and reference positions between one or more :term:`beams <beam>`.

      See also :term:`beam switching`.

    TCAL
      Equivalent temperature of the :term:`noise diode`. Usually given in K.
      This is a function of frequency and time. It is stored as a single value in :term:`SDFITS` as the TCAL column.
      This single value is the average of the noise diode frequency-resolved temperature measurements over the corresponding :term:`IF`.

    VEGAS
      Versatile GBT Astronomical Spectrometer - https://www.gb.nrao.edu/vegas/

    waterfall plot
      A plot (or two-dimensional image) that shows time vs. frequency.

    window
      See :term:`spectral window`.


..    The velocity of a source using the relativistic definition of the velocity-frequency relationship.

..    The velocity of a source using the optical definition of the velocity-frequency relationship.

..    The velocity of a source using the radio definition of the velocity-frequency relationship.


.. _data_org:

Data : Project Code / Session ID
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally projects are assigned a project code, e.g. *AGBT21B-024*, which is
then observed in a number of sessions, numbered starting with 01. The :term:`SDFITS` data associated
with these are stored under ``$SDFITS_DATA``, e.g. for session 5 in this example, this would be
in ``$SDFITS_DATA/AGBT21B_024_05/``.

Possible confusion: a project code "GBT21B-024", is labeled "AGBT21B_024" as the
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
