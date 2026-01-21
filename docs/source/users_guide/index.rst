.. _usersguide:

***********
Users Guide
***********

Overview of dysh
----------------
dysh is a single-dish radio astronomical data reduction package with
emphasis on reducing and analyzing data from the `Green Bank Telescope
<https://greenbankobservatory.org/about/telescopes/gbt>`_.  Our goal in
writing dysh was to create a new, rich Python software environment for
interaction with GBT spectral line data by observers, post-observation
scientists, and GBO staff.

dysh is a toolkit framework designed to not only do standard data
reduction but to put power into users hands by letting them develop
their own sophisticated algorithms, and to make it straightforward for
future developers to add functionality.  The package is designed to allow
fine-grained access to the spectral data and metadata while enabling the
user to also interact with and view the data utilizing simple functions.

We have chosen dysh's user interfaces to minimize barriers to adoption.
These are through the (i)Python shell, Python scripts, Jupyter notebooks
(either through the classic Jupyter Notebook interface or JupyterLabs),
as well as an interactive display.  The `iPython <https://ipython.org/>`_
and `Jupyter Notebook/Lab <https://jupyter.org/>`_ interfaces are powerful,
flexible, and familiar to many in the astronomical community making them
the natural choice over a custom-built scripting language.
Furthermore, the package is built on astropy, numpy, matplotlib, and
pandas with which many astronomers are comfortable.

.. toctree::
   :maxdepth: 2 
   :hidden:
 
   data_structures
   loading_data
   selection
   flagging
   frequencyswitch
   nodding
   positionswitch
   subbeamnod
   calseq
   vane
   iplotter
   custom_baseline
   velocity_frames
   apeff_surf_error
   on_the_fly
   repeated_scans
   using_spectral_weights
   smoothing
   align_spectra
   gauss_fit
   cog
   quality
   line_search
   metadata_management
   dataIO
   merge_sdfits
   hi_survey

