
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
    - GBTFITSLoad.getfs() implemented
    - FSScan implemented

- `ScanBlock <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.scan.ScanBlock>`_  API change
    - `timeaverage() <https://dysh.readthedocs.io/en/release-0.3.0/modules/dysh.spectra.html#dysh.spectra.scan.ScanBlock.timeaverage>`_ now returns a Spectrum instead of a list.  Previously the list contained the time average of each Scan within the ScanBlock.   Now the time average across all Scans in the ScanBlock is returned.

Bug Fixes
-----------
-  `SubBeamNod error when using cycle method <https://github.com/GreenBankObservatory/dysh/issues/207>`_
-  `Spectrum arithmetic operations not working <https://github.com/GreenBankObservatory/dysh/issues/208>`_
-  `SDFITS summary() reports wrong number of integrations <https://github.com/GreenBankObservatory/dysh/issues/211>`_
- `Certain old GBTIDL files could not be read by dysh <https://github.com/GreenBankObservatory/dysh/issues/216>`_
