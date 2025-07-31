.. _glossary:

A dysh glossary
---------------

Commonly used terms.

.. glossary::


    band
      A coherent section of channels in frequency space, all with
      the same channel width. Sometimes called an IF band.
      Explain L-band, K-band, Q-band, W-band etc. Need a table.
      See also **ifnum**

    beam
      The footprint of one receiver horn on the sky. ARGUS has a 
      4x4 multi-beam receiver, numbered 0 through 15.
      Not to be confused with the
      **FWHM**.  At 115 GHz the **FWHM** is about 16", at 86 GHz about
      21".  The beam separation is xx" for ARGUS.
    
      Note that for some instruments beams are also interpreted while
      including other simulteanously taken data in another band/polarization
      See also **fdnum**

    Beam Switching
      This is a variation on position switching using a receiver
      with multiple beams. The "Main" and "Reference" positions on the sky are
      calculated so that the receiver is always pointing at the source. This is most
      useful for point sources.

    blanking vs flagging vs masking
      flagging defines the rules for blanking data, blanking actually applies them

    caloff
      Signal with no calibration in the signal path.

    calon
      Signal with a calibration in the signal path.

    CALSEQ
      Specific W-band calibration procedure involving a sky, cold1 and cold2 measurement.
      See GBT memo 305?
   
    ECSV
      (Enhanced Character Separated Values) a self-describing ascii table format popularized by astropy.
      See also https://github.com/astropy/astropy-APEs/blob/main/APE6.rst

    fdnum
      Feed Number. 0, 1, ...
      Also used as the fdnum= keyword in getXX()
      See also **beam**

    FITS
      (Flexible Image Transport System): the export format
      for data-cube, although there is also a waterfall cube
      (time-freq-pixel) cube available.  Unclear what we will use for
      pure spectra.  **SDFITS** seems overly complex. CLASS needs to
      be supported. 

    flagging vs. blanking

    FWHM
      (Full Width Half Max): the effective resolution of the
      beam if normally given in **FITS** keywords BMAJ,BMIN,BPA.  The
      term **resolution**

    Frequency Switching
      This is a variation on position switching using a receiver
      where the IF is changed. The "Main" and "Reference" positions on the sky are
      calculated so that the receiver is always pointing at the source. This is most
      useful for point sources.    

    horn
      Another term used for :term:`beam` or :term:`pixel`.

    ifnum
      IF number (0,1,...)
      Also used as the ifnum= keyword in getXX()    

    intnum
      Integration number. 0 being the first.
      Also used as the intnum= keyword in getXX()    

    masking vs. flagging blanking

    noise diode
      calibration gadgets at lower frequencies. See also VANE

    OTF Mapping
      On-The-Fly (OTF) mapping is where the telescope is continuesly scanned across the sky
      to sample the emission.
      The samples are then "gridded" into a map, after calibration. Calibration schemes can differ,
      e.g. OnOff in separate scans,  and in-scan definitions.
   
    plnum
      Polarization number (0,1,...). Usually 0 and 1, but of course up to 4 values could be present
      for a full Stokes.
      Also used as the plnum= keyword in getXX()    

    Position Switching
      This is a standard way to obtain spectra by switching
      between a "Main" and "Reference" position on the sky, usually using a single beam. For our
      multi-beam receivers see also Beam Switching

    resolution
      this term is used in the gridder, but it's not
      **FWHM**, it's lambda/D.  Keyword --resolution= is used If
      selected this way, FWHM is then set as 1.15 * resolution. But if
      resolution is chosen larger, what is the effective FWHM?  It
      would be better to have a dimensionless term for
      **resolution/pixel** and a different name for resolution
      alltogether.

    row - often used to refer to a spectrum, technically a row in the fits BINTABLE.

    RRL - Radio Recombination Line

    Ruze's equation
      relating the gain of an antenna to its surface accuracy `eps`
      G = G_0 exp(-(4.pi.eps/lambda)^2)

    Scan - GBT differentiates between different types of scans
     (FSScan, PSScan, TPScan, SubBeamNod Scan). They contain a series of
     calibrated or uncalibrated raw spectra.
   
    ScanBlock - GBT. A container for a series of **scan**'s that looks like a list for the user.
    
    SDFITS
      Single Dish **FITS** format, normally used to store
      raw or even calibrated spectra in a FITS BINTABLE format.  Each
      row in a BINTABLE has an attached RA,DEC (and other meta-data),
      plus the whole spectrum. This standard was drafted in 1995 (Liszt),
      and has been implemented by many telescopes (Arecibo, FAST, GBT, Parkes, ....)

    SFL
      Sanson-Flamsteed projection, used in LMT **FITS** files
      (the GLS - GLobal Sinusoidal is similar to SFL).

    Spectral Window
      In ALMA commonly abbreviated as **spw**, this is closest to what we call a **bank**,
      or **band**, a set of linearly spaced channels.

    Spectrum
      A coherent section in frequency space, with its own unique meta-data (such as polarization,
      ra, dec, time). Normally the smallest portion of data we can assign. A spectrum is
      defined by its own seting of *(crval, crpix, cdelt)* in a FITS WCS sense.
      See also :ref:`storage`.

    SubBeamNod
      Another scan mode

    VANE
      That thing

    VANECAL
      Calibration procedure for ARGUS

    Waterfall Plot
      Tyically

    Window
      See **Spectral Window**

.. _overloaded:

Overloaded Terms
~~~~~~~~~~~~~~~~

Terms used in the code may not exactly match terms used by the develpers of the instruments.
Here we clarify those overloaded terms in the form of a table

.. list-table:: **Table of some overloaded terms**
   :header-rows: 1
   :widths: 15,15,15,45      

   * - code term
     - RSR term
     - SLR term
     - comments
   * - beam
     - pixel?
     - pixel
     - multi-beam receiver
   * - cell
     - n/a
     - cell
     - size of a sky pixel in gridding, usually 2-3 times smaller than the resolution
   * - band
     - board
     - bank
     - spectrometer window
   * - n/a
     - chassis
     - n/a
     - tuple of (pol,beam)
   * - channel
     - channel
     - channel
     - with a simple FREQ WCS{crval,crpix,cdelt}

.. _storage:

Data Dimensions
~~~~~~~~~~~~~~~

This section is not meant to describe either the data format, but the
mental storage model we have in mind to be encapsulated in a Python
class.

A unified data storage of LMT spectra would (naturally) break up the
spectra, such that each spectrum has a different
time, beam, band, polarization, etc.  Each spectrum
can be described as a set of sequential channels, described with a single
*(crval,crpix,cdelt)*) WCS.
In Python row-major array notation where the most slowly varying dimension comes
first this could be written as an **NDarray**:

.. code-block::

      data[ntime, nbeam, npol, nband, nchan]

where we added the ``ntime`` and ``nchan`` as the slowest resp. fastest running dimension
in this row-major (python/C) notation.


.. note:: For those used to GBTIDL **plnum** = **npol**, **ifnum** = **nband**, and
   **fdnum** = **nband**.  Arguably different scans can act as as **ntime**, although
   each scan will often have several snapshots inside of them. ?? **intnum**

.. code-block::

      Overloaded words, including GBT lingo:

      plnum   pol
      fdnum   feed     beam    pixel
      ifnum   window   band

Taking out those an observation can be seen as a set of spectra:

.. code-block::

      spectrum[nbeam, npol, nband]

This exactly matches the concepts used in an SDFITS file, although in the general
definition of SDFITS there is no assumption of the data being able to be stored
in an **NDarray** type array, where the more general

.. code-block::

       sdfits_data[naxis2, ndata]

where in general ``ndata=nchan``, but dialect with ``ndata = npol * nchan`` are
seen in the wild (FAST, Parkes). The FITS name ``naxis2`` is the number of rows,
which is the product of ``time,beam,band,pol`` in our case.


Taking an inventory of current and known future LMT Spectral Line instruments:

* RSR:
  two beams, two pols, 6 bands, though the term *chassis* is used to point at any
  tuple of (beam,pol). So here we have nbeam=2, npol=2,nband=6, nchan=256 and ntime
  it typically 10-20. Each beam happens to look at the same sky position here.

.. note::  If an instrument like RSR would multiplex the (beam,pol) pairs, this would be a challenge
	   to the assumption of homogeneity, and the SDFITS model would be more appropriate.

* 1MM:
  one beam, two pols, two sidebands. So here we have nbeam=1, bpol=2, nband=2, nchan=2k

* SEQ:
  16 beams (though 4 beams per roach board, and each roach board has its own time) in one
  band (they also call it bank) and one polarization. Thus nbeam=16, npol=1, nband=1.
  Once the 2nd IF will be installed, 32 beams will be recognized by the software,
  but organizationally it is easier to to think of 16 beams and 2 bands.

.. note::  The timestamps for the different roach boards make it impossible to store
	   the data in a multi-dimensional array, unless (typicall one) integration
	   is removed. Keeping all data would require ``data[ntime4, 1, 1, 1, nchan]`` for SEQ.

* OMA
  8 beams, 2 bands (banks), 2 polarizations.

* B4R
  4 XFFTS boards, 2.5 GHz/board:  1 beam, 2 bands (USB and LSB), 2 polarizations (XX and YY)

Note that FAST is the only known case that stores data as  ``data[ntime, nchan, npol]``, where
``nchan`` is not the fastest running dimension, but ``npol``. Technically this appears to be the
case such that they can vary ``nchan`` per row.


We thus arrive at the following summary for the multi-dimensional data[] array:

.. code-block::

      data[ntime, nbeam, npol, nband, nchan]

in the table we leave out the ``ntime`` dimension    

.. list-table:: **Table of data dimensions of LMT SLR instruments**
   :header-rows: 1
   :widths: 15,10,10,10,10,30

   * - **data**
     - **nbeam**
     - **npol**
     - **nband**
     - **nchan**
     - comment
   * - RSR
     - 2
     - 2
     - 6
     - 256
     - (pol,beam) tuples are the 4 chassis. 6 overlapping bands make one final spectrum
   * - SEQ
     - 16
     - 1
     - 1 (2)
     - 2k, 4k, 8k
     - beams have time issue, perhaps ntime ~ ntime * nbeam, and nbeam=1. Future will have 2 bands
   * - OMA 
     - 8
     - 2
     - 2
     - 2k, 4k, 8k
     - Future instrument, with 4 more roach boards (USB+LSB)
   * - 1MMRx
     - 1
     - 2
     - 2
     - 2k, 4k, 8k
     - band: 2 IF's in USB/LSB
   * - B4R
     - 1
     - 2
     - 2
     - 32k
     - Japanese 2mm receiver

Single Dish Math
~~~~~~~~~~~~~~~~

The meat of Single Dish math is getting the system temperature


.. math::

   T_{sys} = T_{amb} { { SKY } \over { HOT - SKY } }

or

.. math::

   T_{sys} = T_{cal} { { <SKY> } \over { <HOT - SKY> } } + T_{cal}/2

where the :math:`< >` operator averages over the center 80% of the spectrum.
This way :math:`T_{sys}` is a scalar. The routine ``meantsys`` computes this.

and using this system temperature, calculating the signal by comparing an *ON* and *OFF* position,
assuming there is only sky in the *OFF*:

.. math::

   T_A = T_{sys}  {   { ON - OFF } \over {OFF} }

All of these have values for each channel. How exactly the :math:`T_{sys}` is computed (scalar, vector,
mean/median) is something we generally leave open.

.. math::


is the effective beam (see also GBT memo 296, and gbtpipe/Gridding.py)

.. code-block::

    1.18 * (c / nu0 / 100.0) * 180 / np.pi  # in degrees

Reduction of noise with smoothref=N:

.. math::

     \sigma_N = \sigma_1 \sqrt{   {N+1} \over  {2N}  }

Weight factors

.. math::

      \frac{\Delta t \Delta\nu}{T_{\rm{sys}}^{2}}      \tag{4}

Effective exposure time

.. math::
  { t_{ON} * t_{OFF} } \over {  t_{ON} + t_{OFF}  }

   

Observing: ObsNum / SubObsNum / ScanNum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An observation with a single dish such as LMT is done via proposals, which gets assigned a proposal ID,
associated with the P.I. name. An example of such is **2018-S1-MU-46**

An observation is that divided in a set a **ObsNum** 's, which can be hierchically
divided up in **SubObsNum**'s and **ScanNum**'s. When
an observing script executes, each source will gets its own **ObsNum**, though
calibration data often gets another **ObsNum**.


Band Designations
~~~~~~~~~~~~~~~~~


Nomenclature comes from the IEEE radar band names, but there is also
a NATO nomenclature standard. See e.g. https://en.wikipedia.org/wiki/Radio_spectrum

There is also a list of instruments at GBT clarifying this list
on https://greenbankobservatory.org/portal/gbt/instruments/ or
https://dss.gb.nrao.edu/receivers/summary

* L: 1-2 GHz
* S: 2-4 GHz
* C: 4-8 GHz
* X: 8-12 GHz
* Ku: 12-18 GHz
* K: 18-27 GHz
* Ka: 28-40 GHz
* V: 40-75 (is this used at all)
* Q: 33-50 GHz or 6-9.1mm
* W: 84-116 GHz (also: 75-110)


Suggested Common Parameter Descriptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*   apply_flags : boolean, optional.  If True, apply flags before calibration.


*   bintable : int
            the index for BINTABLE containing the scans

*   channel : number, or array-like
       The channels to operate on
    


*   gbtfits : `~fits.sdfitsload.SDFITSLoad`
        input SDFITSLoad object

*   fdnum : int
            the feed index

*   feeds : int list of two beams
            the two nodding beams

*   fitsindex : int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

*   ifnum: int
            the IF index

*   intnum


*   plnum: int
            the polarization index




*   selection : `~pandas.DataFrame`
            selection object


*   observer_location



*   scans : int or 2-tuple
      The scan(s) to use. A 2-tuple represents (beginning, ending) scans.
      Default: show all scans

*   scans : array-like
            list of one or more scans
   
   

*   overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.
*    scan: int        # TPScan(
        scan number


*   scans : dict        # PSScan
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference)
        scan numbers NOTE: there should be one ON and one OFF, a pair

*   scan : dict     # NodScan
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference) scan numbers
        NOTE: there should be one ON and one OFF, a pair. There should be at least two beams (the nodding beams)
        which will be resp. on source in each scan.

*   scan : int    # FSScan(ScanBase)
        Scan number that contains integrations with a series of sig/ref and calon/caloff states.

*   smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration


*   weights: str
            'tsys' or None.
    

Finally
~~~~~~~

No more.

   
