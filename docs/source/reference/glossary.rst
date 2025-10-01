.. _glossary:

A dysh glossary
---------------


.. glossary::


    argus
      yes, this is an abbreviation of sorts - https://www.gb.nrao.edu/argus/

    band
      A coherent section of channels in frequency space, all with
      the same channel width. Sometimes called an IF band.
      See also :term:`ifnum`

    bank
      Overloaded term for a **band**, possibly referring to hardware.

    beam
      The footprint of one receiver horn on the sky. ARGUS has a
      4x4 multi-beam receiver, numbered 0 through 15.
      Not to be confused with the
      **FWHM**.  At 115 GHz the **FWHM** is about xx", at 86 GHz about
      xx".  The beam separation is xx" for ARGUS.

      Note that for some instruments beams are also interpreted while
      including other simulteanously taken data in another band/polarization
      See also :term:`fdnum``

      See also :term:`multi-beam`

    Beam Switching
      This is a variation on position switching using a receiver
      with multiple beams. The "Main" and "Reference" positions on the sky are
      calculated so that the receiver is always pointing at the source. This is most
      useful for point sources. See also :term:`Position Switching`

    blanking vs flagging vs masking
      flagging defines the rules for blanking data, blanking actually applies them

    caloff
      Signal with no calibration in the signal path.

    calon
      Signal with a calibration in the signal path.

    DYSH_DATA
      (optional) environment variable pointing to a directory for a convenient view of
      data for developers.
      See also :term:`SDFITS_DATA`.

    ECSV
      (Enhanced Character Separated Values) a self-describing ascii table format popularized by astropy.
      See also https://github.com/astropy/astropy-APEs/blob/main/APE6.rst

    fdnum
      Feed Number. 0, 1, ...
      Also used as the **fdnum=** keyword in getXX()
      See also :term:`beam`

    FITS
      (Flexible Image Transport System): the export format
      for data-cube, although there is also a waterfall cube
      (time-freq-pixel) cube available.  Unclear what we will use for
      pure spectra.  **SDFITS** seems overly complex. CLASS needs to
      be supported.

    flagging vs. masking
      flagging is non-destructive.

    FWHM
      (Full Width Half Max): the effective resolution of the
      beam if normally given in **FITS** keywords BMAJ,BMIN,BPA.  The
      term **resolution**

    Frequency Switching
      This is a variation on position switching using a receiver
      where the IF is changed. The "Main" and "Reference" positions on the sky are
      calculated so that the receiver is always pointing at the source. This is most
      useful for point sources.

    getXX()
      Generic name for the dysh access routines, e.g. getps, getfs, getnod etc.

    horn
      Another term used for :term:`beam` or :term:`pixel`.

    ifnum
      IF number (0,1,...)
      Also used as the ifnum= keyword in getXX().
      See also :term:`band` and :term:`window`

    intnum
      Integration number. 0 being the first.
      Also used as the intnum= keyword in getXX()

    kfpa
      K-band Focal Plane Array

    masking vs. flagging vs. blanking
      Masking removes or hides pixels,
      while flagging keeps the pixels but attaches a status to them for later filtering or analysis. (google)

      blanking is destructive.

      OK  google is also very conflicted here.  Compare python:

      In python a mask is True/False, where True indicates an element of the array is to be selected.

    multi-beam
      If an instrument has multiple beams that typically point are different areas in the sky
      (e.g. **ARGUS** in a 4x4 configuration, and **Kfpa** in a 7 beam hexagonal shape).

    Nod or Nodding
      An observing mode ...

    OTF Mapping
      In this procedure the telescope is scanned across the sky to sample the emission.
      The samples are then "gridded" into a map.

    pixel
      An overloaded term. Sometimes referred to as the :term:`beam`, but usually interpreted as
      the size of a single (usually square) element in a gridded map (e.g. from an OTF), which
      we commonly also refer to as a *picture element*.

    plnum
      Polarization number (0,1,...). Usually 0 and 1, but of course up to 4 values could be present
      for a full Stokes.
      Also used as the plnum= keyword in getXX()

    Position Switching
      This is a standard way to obtain spectra by switching
      between a "Main" and "Reference" position on the sky, usually using a single beam. For our
      multi-beam receivers see also Beam Switching

    Project ID
      Or what's the name at GBO?

    resolution
      this term is used in the gridder, but it's not
      **FWHM**, it's lambda/D.  Keyword --resolution= is used If
      selected this way, FWHM is then set as 1.15 * resolution. But if
      resolution is chosen larger, what is the effective FWHM?  It
      would be better to have a dimensionless term for
      **resolution/pixel** and a different name for resolution
      alltogether.

    RRL - Radio Recombination Line
      more tbd

    Scan - GBT differentiates between different types of scans
     (FSScan, PSScan, TPScan, SubBeamNod Scan). Each of these comes
     with a corresponding :term:`getXX()`

    ScanBlock
      A container for a series of **scan**'s.
      There is a whole section on this. explanations/scanblock/index.html

    SDFITS
      Single Dish **FITS** format, normally used to store
      raw or even calibrated spectra in a FITS BINTABLE format.  Each
      row in a BINTABLE has an attached RA,DEC (and other meta-data),
      plus the whole spectrum. This standard was drafted in 1995 (Liszt),
      and has been implemented by many telescopes (Arecibo, FAST, GBT, Parkes, ....)

      There are two other links in the manual: reference/sdfits_files/index.html
      and   explanations/sdfits/index.html

    SDFITS_DATA
      (optional) environment variable pointing to a directory where SDFITS files
      and projects are stored.

    SFL
      Sanson-Flamsteed projection, sometimes used in gridding OTF maps.
      (the GLS - GLobal Sinusoidal is similar to SFL).

    Spectral Window
      In ALMA commonly abbreviated as **spw**, this is closest to what we call a **bank**,
      or **band**, a set of linearly spaced channels. See also :term:`ifnum`

    Spectrum
      A coherent section in frequency space, with its own unique meta-data (such as polarization,
      ra, dec, time). Normally the smallest portion of data we can assign. A spectrum is
      defined by its own seting of *(crval, crpix, cdelt)* in a FITS WCS sense.

    SubBeamNod
      Subreflect Beam Nodding. The getXX() is now called `subbeamnod`

    VEGAS
      Versatile GBT Astronomical Spectrometer - https://www.gb.nrao.edu/vegas/

    Window
      See **Spectral Window**


Band Designations
~~~~~~~~~~~~~~~~~


W-band

Q-band

L-band

K-band


.. _sdmath:



Single Dish Math
~~~~~~~~~~~~~~~~

The meat of Single Dish math is getting the system temperature


.. math::  :label: sdmath1

   T_{sys} = T_{amb} { { SKY } \over { HOT - SKY } }

or

.. math:: :label: sdmath2

   T_{sys} = T_{cal} { { <SKY> } \over { <HOT - SKY> } } + T_{cal}/2

where the :math:`< >` operator averages over the center 80% of the spectrum.
This way :math:`T_{sys}` is a scalar. The routine ``meantsys`` computes this.

and using this system temperature, calculating the signal by comparing an *ON* and *OFF* position,
assuming there is only sky in the *OFF*:

.. math:: :label: sdmath3

   T_A = T_{sys}  {   { ON - OFF } \over {OFF} }

All of these have values for each channel. How exactly the :math:`T_{sys}` is computed (scalar, vector,
mean/median) is something we generally leave open.

.. math::


the effective beam (see also GBT memo 296, and gbtpipe/Gridding.py - at 109 GHz)

.. code-block::

    1.18 * (c / nu0 / 100.0) * 180 / np.pi  # in degrees

Reduction of noise with smoothref=N:

.. math:: :label: sdmath4

     \sigma_N = \sigma_1 \sqrt{   {N+1} \over  {2N}  }

Weight factors are 1/:math:`\sqrt(RMS)`

.. math:: :label: sdmath4

      \frac{\Delta t \Delta\nu}{T_{\rm{sys}}^{2}}

Effective exposure time

.. math:: :label: sdmath5

  { t_{ON} * t_{OFF} } \over {  t_{ON} + t_{OFF}  }

As shown in :eq:`sdmath2` we can ...


Data : Project ID / Session / Scan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally projects are assigned a project id, e.g. *AGBT21B_024*, which is
then observed in a number of sessions, numbered starting with 1. The SDFITS data associated
with these are stored under **$SDFITS_DATA**, e.g. for session 5 of the example above, this would be
in **$SDFITS_DATA/AGBT21B_024_05/**.   At GBO  SDFITS_DATA=/home/sdfits, but outside
of GBO this will be user defined. Another default is **$DYSH_DATA/sdfits**, if
**DYSH_DATA** is used.
