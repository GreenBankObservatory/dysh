**********************
dysh Documentation
**********************

dysh is a Python spectral line data reduction and analysis program for single dish data with specific emphasis on data from the Green Bank Telescope.
It is currently under development in collaboration between the `Green Bank Observatory <https:/greenbankobservatory.org>`_ and the Laboratory for Millimeter-Wave Astronomy (LMA) at the `University of Maryland (UMD) <https://www.astro.umd.edu>`_.
It is intended to replace `GBTIDL <https://gbtidl.nrao.edu/>`_, GBO's current spectral line data reduction package.

Installation
============

.. code:: bash

    pip install "dysh[nb] @ git+https://github.com/GreenBankObservatory/dysh"

Read on for a brief overview.

.. tip::

    You can find detailed installation instructions in :ref:`installing-dysh`.

Launching dysh
==============

After being installed, the ``dysh`` command will be available through the command line.
This will launch an iPython session with some modules and classes pre-loaded (e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad`), and with logging.
We refer to this interface as the dysh shell.

Loading Data
============

dysh can read and write :ref:`SDFITS files <sdfits-explanation>`.
Once inside the dysh shell loading an :ref:`SDFITS file <sdfits-explanation>` can be done using:

.. code:: Python

    from dysh.util.download import from_url
    # The following line will download 31 MB to the current working directory.
    filename = from_url("http://www.gb.nrao.edu/dysh/example_data/frequencyswitch/data/TREG_050627/TREG_050627.raw.acs/TREG_050627.raw.acs.fits")
    sdfits = GBTFITSLoad(filename) # This will load the SDFITS file(s) at `filename`.

In the above code, you can replace ``filename`` with a path to your own data.
Either a single SDFITS file or a directory containing multiple SDFITS files.

The contents of the loaded data can be inspected using the `GBTFITSLoad.summary <dysh.fits.gbtfitsload.GBTFITSLoad.summary>` method, as:

.. code:: Python

    sdfits.summary()

.. raw:: html
    :file: files/example_getfs_summary.html

This shows five scans of W3OH.
We know that this data was observed using :term:`frequency switching`.

Calibrating Data
================

Since the data was observed using frequency switching, we use the `GBTFITSLoad.getfs <dysh.fits.gbtfitsload.GBTFITSLoad.getfs>` method to calibrate the data.
This method requires the :term:`fdnum`, :term:`ifnum` and :term:`plnum` parameters, and it will return a :ref:`ScanBlock <scanblocks>` object.

.. code:: Python

    scan_block = sdfits.getfs(fdnum=0, ifnum=0, plnum=0)

That's it, now ``scan_block`` contains all the calibrated data.


The calibrated data can be plotted as a waterfall using `ScanBlock.plot <dysh.spectra.scan.ScanBlock.plot>`.
To plot the data as a waterfall saturating the color scale at 10 K use:

.. code:: Python

    scan_block_plot = scan_block.plot(vmax=10)

.. image:: files/scan_block_plot_vmax10.png


This method returns a `~dysh.plot.scanplot.ScanPlot` object, which can be used to modify the plot.
See the :ref:`waterfall plot guide <scanblock-plots>` for more details.

Time Averaging Data
===================

Time averaging the calibrated data into a single spectrum can be done using the `ScanBlock.timeaverage <dysh.spectra.scan.ScanBlock.timeaverage>` method, like:

.. code:: Python

    spectrum = scan_block.timeaverage()

The return from this method is a `~dysh.spectra.spectrum.Spectrum` object.

Plot the spectrum using the `Spectrum.plot <dysh.spectra.spectrum.Spectrum.plot>` method:

.. code:: Python

    spec_plot = spectrum.plot()

.. image:: files/spectrum_plot.png


This method returns a `~dysh.plot.specplot.SpectrumPlot` object, which can be used to modify the plot.
See the :ref:`spectrum plot guide <spectrum-plots>` for more details.

Baseline Removal
================

To remove a baseline from the data use the `Spectrum.baseline <dysh.spectra.spectrum.Spectrum.baseline>` method.
This can use different polynomials to fit and remove a baseline from the spectrum.
To fit and remove a Chebyshev polynomial of degree 17, ignoring regions with spectral lines, and their ghosts, use:

.. code:: Python

    # Define the channel ranges to be
    # excluded from the baseline fit.
    exclude = [(12720, 13259),
               (14060, 14440),
               (16000, 16420),
               (17375, 17640),
               (19355, 19660),
               (20575, 20955)]
    spectrum.baseline(model="chebyshev", degree=17, remove=True, exclude=exclude)

The spectrum has been baseline removed.
Plot again, but focus on the line-free regions:

.. code:: Python

    spec_plot = spectrum.plot(ymax=1, ymin=-1)

.. image:: files/spectrum_plot_bsub.png

Slicing a Spectrum
==================

The brightest line is the OH maser at 1665.4018 MHz.
To crop the spectrum around this line we use a `slice`.
We will use `astropy.units` to define the frequency range:

.. code:: Python

    from astropy import units as u
    oh_bright_spec = spectrum[1665.3*u.GHz:1665.9*u.MHz]

Change the rest frequency to that of the 1665.4018 MHz line:

.. code:: Python

    oh_bright_spec.doppler_rest = 1665.4018 * u.MHz

Plot in velocity units (``xaxis_unit``) using the local standard of rest as velocity frame (``vel_frame``) and the radio Doppler convention (``doppler_convention``):

.. code:: Python

    oh_bright_spec.plot(xaxis_unit="km/s", vel_frame="lsrk", doppler_convention="radio")

.. image:: files/oh_bright_kms_lsrk_radio.png

Note that the setting the units, velocity frame or Doppler convention during plotting does not modify the spectrum.

Saving a Spectrum
=================

Save the spectrum to a FITS file:

.. code:: Python

    oh_bright_spec.write("W3OH_1665MHz_pol0.fits", format="fits")


Reporting Issues
================

If you find a bug or something you think is in error, please report it on
the `GitHub issue tracker <https://github.com/GreenBankObservatory/dysh/issues>`_.
You must have a `GitHub account <https://github.com>`_ to submit an issue.


Contents
===============

.. grid:: 1 2 2 2

    .. grid-item-card::
        :shadow: md
        :margin: 2 2 0 0

        :octicon:`mortar-board;3em;orange` **Tutorials**

        Learning-oriented lessons take you through a series
        of steps to complete a project.

        Most useful when you want to get started reducing your data.

        .. button-link:: tutorials/index.html
            :color: primary
            :outline:
            :click-parent:

            Go to Tutorials

    .. grid-item-card::
        :shadow: md
        :margin: 2 2 0 0

        :octicon:`terminal;3em;green` **Recipes**

        Practical step-by-step guides to help you achieve a specific goal.

        Most useful when you're trying to get something done.


        .. button-link:: how-tos/index.html
            :color: primary
            :outline:
            :click-parent:

            Go to Recipes

    .. grid-item-card::
        :shadow: md
        :margin: 2 2 0 0

        :octicon:`repo;3em;purple` **Explanation**

        Big-picture explanations of higher-level concepts.

        Most useful for building understanding of a particular topic.


        .. button-link:: explanations/index.html
            :color: primary
            :outline:
            :click-parent:

            Go to Explanation Material

    .. grid-item-card::
        :shadow: md
        :margin: 2 2 0 0

        :octicon:`tools;3em;sd-text-primary` **References**

        Nitty-gritty technical descriptions of how `dysh` works.

        Most useful when you need detailed information about the API or how to
        contribute.


        .. button-link:: reference/index.html
            :color: primary
            :outline:
            :click-parent:

            Go to Reference Material


.. toctree::
   :maxdepth: 2
   :hidden:

   whatsnew/CHANGES.rst
   getting_started/index
   tutorials/index
   how-tos/index
   explanations/index
   reference/index
   for_beta_testers/index
   for_developers/index
   reference/glossary


Credits
=======
dysh is being developed by a partnership between the Green Bank Observatory and the Laboratory for Millimeter-wave Astronomy at the University of Maryland, College Park.

Dev Team
--------
| Marc Pound (UMD)
| Peter Teuben (UMD)
| Pedro Salas (GBO)
| Evan Smith (GBO)
| Thomas Chamberlin (GBO)

Earlier contributors
--------------------
| Victoria Catlett (former GBO)
