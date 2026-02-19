
.. _usersguide-launching:

**************
Starting dysh
**************

dysh shell
==========

dysh comes with two pre-configured interfaces that can be started from the
command line. The first, ``dysh``, will start an `IPython <https://ipython.readthedocs.io/>`_ shell with common
dysh classes, and astropy, numpy, and pandas modules already imported.
``dysh`` has a number of command line options, viewable with ``--help``

.. _usersguide-launching-dysh:

.. code:: bash

   dysh --help

   usage: dysh [-h] [-i] [-p PATHS [PATHS ...]] [-P PROFILE] [-L FITS_LOADER] [--colors {NoColor,Neutral,Linux,LightBG}] [-v {0,1,2,3}]
            [--log LOG] [-q] [--version] [--skip-config] [--hide-tb]
            [file]

   Dysh interactive shell. All CLI arguments other than those defined below are passed through to ipython; see $ ipython --help for more
   details

   positional arguments:
     file                  Path to a Dysh Script to run

   options:
     -h, --help            show this help message and exit
     -i, --interactive     Remain in interactive mode after running a script
     -p PATHS [PATHS ...], --paths PATHS [PATHS ...]
                           FITS file paths to load initially
     -P PROFILE, --profile PROFILE
                           The IPython profile to use
     -L FITS_LOADER, --fits-loader FITS_LOADER
                           The SDFITS loader class name to use
     --colors {NoColor,Neutral,Linux,LightBG}
                           Set the color scheme
     -v {0,1,2,3}, --verbosity {0,1,2,3}
                           Set logging verbosity
     --log LOG             Specify log path
     -q, --quiet           Silence DEBUG- and INFO-level logs to stderr
     --version             Print version and exit
     --skip-config         Skip creating a configuration file
     --hide-tb             Hide traceback

.. _usersguide-launching-dysh-lab:

dysh lab
========

There is also a custom `JupyterLab <https://jupyterlab.readthedocs.io/>`_ interface, ``dysh-lab``, which will start a JupyterLab server and open a launcher in your browser.

.. note::
    At GBO ``dysh-lab`` will not launch a web browser. You must direct your web browser to the appropriate location, and `set up a tunnel <https://gbtdocs.readthedocs.io/en/latest/how-tos/infrastructure/remote-connection.html#vnc-connection>`_ if working through ssh.

.. code:: bash

   dysh-lab --help

   usage: dysh-lab [-h] [--version] [--no-browser]

   Dysh lab

   options:
     -h, --help    show this help message and exit
     --version     Print version and exit
     --no-browser  Do not open browser automatically
