
Version 0.4.0
=============

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

Functionality & Enhancements
---------------------------------
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

Bug Fixes
-----------
-  `SubBeamNod error when using cycle method <https://github.com/GreenBankObservatory/dysh/issues/207>`_
-  `Spectrum arithmetic operations not working <https://github.com/GreenBankObservatory/dysh/issues/208>`_
-  `SDFITS summary() reports wrong number of integrations <https://github.com/GreenBankObservatory/dysh/issues/211>`_
- `Certain old GBTIDL files could not be read by dysh <https://github.com/GreenBankObservatory/dysh/issues/216>`_
