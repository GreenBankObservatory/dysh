# Example Notebooks

This directory contains Jupyter notebooks with ``dysh``
examples. These (and more) are all compiled into `readthedocs`, so we
normally recommmend reading the ``dysh`` documentation at
https://dysh.readthedocs.io , which contains the same contents
organized and with more context.

However, if you want to skip directly into the notebooks, start with `AAAREADME.ipynb` as a top level guide

## Notebook layout

Each notebook has a specific layout that should be followed:

1. Title with small description what the notebook covers
2. Optional Background (nodding has a good example)
3. Dysh Commands gives summary of commands (re)introduces here. Generally function arguments can be skipped, but the
   `velocity_frames`  notebook has an counter example
3. Loading Modules loaded and setup for logging/file I/O
4. Data Retrievel, data should be obtained via 	`dysh_data()`
5. 
9. Final Stats at the end should take a spectrum and the `stats()` -  RMS be recorded so it can be compared
   visually during runtime.
9. See Also - references to related and relevant notebooks. Optional.   

Also:

1. Embedded links should be on start of a new line (they are long, easier to find and change this way). E.g.

    [Spectrum](https://dysh.readthedocs.io/en/latest/modules/dysh.spectra.html#module-dysh.spectra.spectrum)

2. If the last statement in a cell is a `plot()`, it should become `plot();` to prevent an ugly looking
   and useless return object printed to the file.




## Notebook Description

| Notebook.ipynb       | Description                                                               |
|----------------------|---------------------------------------------------------------------------|
| AAAREADME            | Top level guidance in which order to read these notebooks                 |
| align_spectra        | How to align spectra                                                      |
| apeff_surferr        | How to modify aperture efficiency or surface error when calibrating       |
| calseq               | Calibration of GBT W-Band observations                                    |
| custom\_baseline     | How to fit a custom baseline model                                        |
| dataIO               | How to read and save data                                                 |
| flagging             | Data flagging                                                             |
| frequencyswitch      | Calibration of frequency switched observations                            |
| gauss\_fit           | Gaussian Fitting                                                          |
| hi\_survey           | Calibration of L-Band observations using position switching               |
| line\_search         | How to search for potential spectral lines                                |
| merge\_sdfits        | Merging SDFITS Files                                                      |
| metadata\_management | How to access and modify metadata for an SDFITS                           |
| nodding              | Calibration of nod observations                                           |
| on\_the\_fly         | Calibration of on-the-fly observations using position switching at L-Band |
| positionswitch       | Calibration of position switched observations                             |
| quality              | Data Quality                                                              |
| repeated\_scans      | Working With Repeated Scan Numbers                                        |
| selection            | How to select data                                                        |
| smoothing            | How to smooth spectra                                                     |
| subbeamnod           | Calibration of subreflector beam nodding observations                     |
| using_spectral_weights | How to set and use weights in spectral averaging                        |
| vane                 | Calibration of Argus observations with a vane                             |
| velocity_frames      | How to manipulate spectral axes                                           |

## notebooks alphabetical

This also keeps a log of what's not perfect and needs to be worked on

```
align_spectra
apeff_surferr
calseq                          radiometer is bad
custom_baseline
dataIO

flagging
frequencyswitch
gauss_fit
hi_survey
line_search

merge_sdfits
metadata_management
nodding                          logger,   radiometer is bad
on_the_fly
positionswitch

quality
repeated_scans
selection
smoothing 
staff_training                          --  KeyError: 'TSCALE'

subbeamnod
using_spectral_weights
vane                                  -- W Warning: 0 != 4: inconsistency counters in mask usage
velocity_frames
```
