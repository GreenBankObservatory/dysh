Version 0.7.0
=============

.. _0.7.0-functionality:

Functionality & Enhancements
----------------------------
- Scan/ScanBlock baseline subtraction
  - It is now possible to subtract a baseline model from all the integrations in a `Scan` or `ScanBlock` using the `subtract_baseline` method.
- Curve of Growth
  - `Spectrum` objects now have a `cog` method that computes the curve of growth.

.. _0.7.0-bugfixes:

Bug Fixes
---------
- Fixes a bug where using `proc=...` with `getps` would result in an error.
- Fixes a bug where unpaired OnOff or OffOn observing procedures would result in an error.
- Fixes a bug where using an exclude region with a boundary equal to the number of channels would produce an error.

Version 0.6.2
=============

.. _0.6.2-bugfixes:

Bug Fixes
---------
- Fixes a bug where `Spectrum.average` would not mask the data, if there was anything to be masked.
- Fixes a bug where `Spectrum.average` would not propagate the data reduction history to the resulting `Spectrum`.
- Fixes a bug where the integration number for files with multiple binary tables would have bogus values. See `Issue #425 <https://github.com/GreenBankObservatory/dysh/issues/425>`_.

Version 0.6.1
=============

.. _0.6.1-bugfixes:

Bug Fixes
---------
- Fixes a bug where calibration routines would not fetch the correct data when the input SDFITS had multiple binary tables.
- Fixes a bug where calibration routines would not work if there was no noise diode being fired. This is still an issue for `getfs`.

Version 0.6.0
=============

.. _0.6.0-functionality:

Functionality & Enhancements
----------------------------
- Baseline plotting
  - `Spectrum.plot` will now update when using `Spectrum.baseline`. If `remove=False` the plot will show the baseline solution, if `remove=True` the plot will update to show the baseline subtracted data.

Version 0.5.0
=============

.. _0.5.0-functionality:

Functionality & Enhancements
----------------------------
- Calibration with a flexible Off
  - Implements `GBTFITSLoad.getsigref <https://dysh.readthedocs.io/en/latest/reference/modules/dysh.fits.html#dysh.fits.gbtfitsload.GBTFITSLoad.getsigref>`_.
- Calibration of W-Band and Argus observations
  - There are now `GBTFITSLoad.calseq` and `GBTFITSLoad.vanecal` functions to derive system temperatures for W-Band and Argus observations, respectively.

Version 0.4.0
=============

.. _0.4.0-functionality:

Functionality & Enhancements
----------------------------
- Calibration of Nod observations
    - Nod observations can be calibrated using `GBTFITSLoad.getnod() <https://dysh.readthedocs.io/en/release-0.4.0/reference/modules/dysh.fits.html#dysh.fits.gbtfitsload.GBTFITSLoad.getnod>`_. By default `dysh` will identify which beams were used for the nodding, but these can be specified via the `fdnum` argument.
- Metadata management
    - `SDFITSLoad` and `GBTFITSLoad` objects now have access to their column data via their get and set methods. For an example see the `metadata management recipes <https://dysh.readthedocs.io/en/release-0.4.0/how-tos/examples/metadata_management.html>`_.
- Smoothing
    - `Spectrum` objects can now be smoothed using their `smooth <https://dysh.readthedocs.io/en/latest/release-0.4.0/modules/dysh.spectra.html#dysh.spectra.spectrum.Spectrum.smooth>`_ method.
    - Current smoothing kernels include a Gaussian, boxcar and Hanning windows. Their widths are specified in channels.
    - By default the `Spectrum` will be decimated by the kernel width. This can be changed with the `decimate` argument.
- Alignment of `Spectrum` objects
    - `Spectrum` objects can now be aligned with respect to each other to match their spectral axes.
    - The alignment matches the first channel of the `Spectrum` being aligned.
    - The alignment can be done in frequency or velocity space.
- Data IO
    - `dysh` can now write and read `SDFITSLoad`, `ScanBlock` and `Spectrum` objects to a variety of formats. For an example see `data IO recipes <https://dysh.readthedocs.io/en/release-0.4.0/how-tos/examples/dataIO.html>`_.
- Logging
    - `dysh` now uses `logging`. The logging level can be set at startup using the `-v` option or at runtime using `dysh.log.init_logging <https://dysh.readthedocs.io/en/release-0.4.0/reference/modules/dysh.log.html#dysh.log.init_logging>`_.
- Data reduction history
    - `dysh` objects now keep track of their history in their `history` attribute. This is written to disk to enhance data reduction reproducibility.
- Online/Offline mode
    - `dysh` can now access the latest spectral line observations using `GBTOnline <https://dysh.readthedocs.io/en/release-0.4.0/reference/modules/dysh.fits.html#dysh.fits.gbtfitsload.GBTOnline>`_. This will also automatically update the contents of the `GBTFITSLoad` object as new data is written to disk. This mode is only available when working from a `GBO data reduction host <https://greenbankobservatory.org/portal/gbt/processing/#data-reduction-machines>`_.
    - `GBTOffline <https://dysh.readthedocs.io/en/release-0.4.0/reference/modules/dysh.fits.html#dysh.fits.gbtfitsload.GBTOffline>`_ will fetch the data for an observations by specifying the session id (e.g., AGBT24A_999_01).
- Drop support for `Python3.8` and `Python3.9`
    - `dysh` now is only available for `Python3.10` and above.
- Gain correction classes
    - Gain correction classes for computing gain correction as a function of elevation, aperture efficiency, surface error, and airmass. (See `GBTGainCorrection <https://dysh.readthedocs.io/en/release-0.4.0/reference/modules/dysh.util.html#dysh.util.gaincorrection.GBTGainCorrection>`_).

.. _0.4.0-bugfixes:

Bug Fixes
---------
- `Velocity frame documentation <https://github.com/GreenBankObservatory/dysh/issues/303>`_
- `Error with frame conversion of Spectrum <https://github.com/GreenBankObservatory/dysh/issues/401>`_
- `pytest does not remove the files it creates <https://github.com/GreenBankObservatory/dysh/issues/369>`_
- `specutils excise_regions are not inclusive on the first boundary <https://github.com/GreenBankObservatory/dysh/issues/378>`_
- `use of Table.loc fails for astropy 6.1.0 <https://github.com/GreenBankObservatory/dysh/issues/245>`_
- `Descriptive error/warning message for blank integrations <https://github.com/GreenBankObservatory/dysh/issues/254>`_
- `Notebook download outputs have white background in dark mode <https://github.com/GreenBankObservatory/dysh/issues/336>`_
- `ReadTheDocs raises new warnings <https://github.com/GreenBankObservatory/dysh/issues/338>`_
- `Use of Pathlib <https://github.com/GreenBankObservatory/dysh/issues/347>`_
- `Plotting changes spectral_axis of a Spectrum <https://github.com/GreenBankObservatory/dysh/issues/372>`_
- `gettp() does not separate the IF's in a nodding example <https://github.com/GreenBankObservatory/dysh/issues/361>`_
- `Spectrum smooth does not preserve vel frame <https://github.com/GreenBankObservatory/dysh/issues/417>`_

Version 0.3.0
==============

.. _0.3.0-functionality:

Functionality & Enhancements
----------------------------
- Handling of Doppler frames and conventions
    - A GBT Spectrum has a spectral axis derived from the WCS of the spectrum (which in turn is created from meta data in the SDFITS file).  The default Doppler frame for this axis is topocentric.  The spectral axis of a Spectrum can be converted to standard frames recognized by astropy: LSRK, HCRS, ICRS, GCRS, ITRS, GalactoCentric. See `Spectrum.set_frame <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.spectrum.Spectrum.set_frame>`_, `Spectrum.with_frame <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.spectrum.Spectrum.with_frame>`_, and also the `xaxis_unit` and `vel_frame` keywords to `SpecPlot.plot.  <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.plot.html#dysh.plot.specplot.SpectrumPlot.plot>`_
    -  The Doppler conventions *radio, optical, relativistic* are recognized by `dysh`.  Users can convert a Spectrum to different conventions with `Spectrum.set_convention <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.spectrum.Spectrum.set_convention>`_ and `Spectrum.with_velocity_convention <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.spectrum.Spectrum.set_convention>`_.  See also, the `doppler_convention` keyword of  `SpecPlot.plot <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.plot.html#dysh.plot.specplot.SpectrumPlot.plot>`_
- Data Selection
    - The `Selection <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.util.html#dysh.util.selection.Selection>`_ class implements a very flexible way of selecting data rows from an SDFITS file using any column name.  (Column name aliases are also supported).  Multiple selection rules are logically combined to a final selection.
    - Data selection is implemented on `GBTFITSLoad <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.fits.html#module-dysh.fits.gbtfitsload>`_ via delegation to a Selection attribute.
- Frequency switching calibration
    - FS data can be calibrated using `GBTFITSLoad.getfs() <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.fits.html#dysh.fits.gbtfitsload.GBTFITSLoad.getfs>`_  with the option to fold the signal and reference spectra.
    - The Scan class for calibrating frequency switching, `FSScan <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.scan.FSScan>`_, has been implemented.  Users should not need to create these directly, but rather through *getfs()*.

- `ScanBlock <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.scan.ScanBlock>`_  API change
    - `timeaverage() <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.scan.ScanBlock.timeaverage>`_ now returns a Spectrum instead of a list.  Previously the list contained the time average of each Scan within the ScanBlock.   Now the time average across all Scans in the ScanBlock is returned.

.. _0.3.0-bugfixes:

Bug Fixes
---------
-  `SubBeamNod error when using cycle method <https://github.com/GreenBankObservatory/dysh/issues/207>`_
-  `Spectrum arithmetic operations not working <https://github.com/GreenBankObservatory/dysh/issues/208>`_
-  `SDFITS summary() reports wrong number of integrations <https://github.com/GreenBankObservatory/dysh/issues/211>`_
- `Certain old GBTIDL files could not be read by dysh <https://github.com/GreenBankObservatory/dysh/issues/216>`_
