*******************
Installing ``dysh``
*******************

``dysh`` requires Python 3.8+ and recent versions of
`astropy <https://astropy.org>`_, 
`numpy <https://numpy.org>`_,
`scipy <https://scipy.org>`_, 
`pandas <https://pandas.pydata.org>`_,
`specutils <https://specutils.readthedocs.io/en/stable/>`_,  and
`matplotlib <https://matplotlib.org>`_.

With `pip` from PyPi
====================

``dysh`` is most easily installed with ``pip``, which will take care of
any dependencies.  The packaged code is hosted at the `Python Packaging
Index <https://pypi.org/project/dysh>`_.

.. code::

    $ pip install dysh

From github
===========

To install from github without creating a separate virtual environment: 

.. code::

    $ git clone git@github.com:GreenBankObservatory/dysh.git
    $ cd dysh
    $ pip install -e .

If you wish to install using a virtual environment, which we strongly recommend if you plan to contribute to the code, see `Development`_.

Development
===========

Here are the steps if you want to develop code for ``dysh``. 

Installation
------------

We use `hatch <https://hatch.pypa.io/>`_ to manage the build environment.
The usual caveats apply how you set up your python development environment.

#.  Clone the repo and install hatch.

    .. code-block:: bash

        $ git clone git@github.com:GreenBankObservatory/dysh.git
        $ cd dysh
        $ pip install hatch


#.  Create and activate a virtual environment with hatch and install the packages required for development.
The virtual environment will be created the first time; subsequent invoking ``hatch shell`` will simply load the created environment.

    .. code-block:: bash

        $ hatch shell
        (dysh) $ pip install -r requirements_dev.txt


#.  Build and install the package

    .. code-block:: bash

        (dysh) $ hatch build
        (dysh) $ pip install -e .

#.  You can exit this environment (which effectively had started a new shell).

    .. code-block:: bash

        (dysh) $ exit
        $ 

#.  Each time when you come back in this directory without being in this virtual environment, you'll need to load the virtual environment.

    .. code-block:: bash

        $ hatch shell

Notice you can *only* do that from within the original install directory tree.

Testing
-------
 We use pytest for unit and integration testing.  From the top-level dysh directory, run:

.. code-block:: bash

    $ pytest

