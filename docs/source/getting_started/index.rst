***************
Getting Started
***************

Installing ``dysh``
===================

``dysh`` requires Python 3.9+ and recent versions of
`astropy <https://astropy.org>`_,
`numpy <https://numpy.org>`_,
`scipy <https://scipy.org>`_,
`pandas <https://pandas.pydata.org>`_,
`specutils <https://specutils.readthedocs.io/en/stable/>`_,  and
`matplotlib <https://matplotlib.org>`_.

We strongly recommend the use of a virtual environment for installing `dysh`.

With `pip` from PyPi
--------------------

``dysh`` is most easily installed with ``pip``, which will take care of
any dependencies.  The packaged code is hosted at the `Python Packaging
Index <https://pypi.org/project/dysh>`_.

.. code::

    $ pip install dysh

.. warning::
    `dysh` is currently in development and the above command will install the latest stable version of `dysh` which might not reflect the contents of the documentation.
    For beta testing please see :ref:`beta-install`.

From GitHub
-----------

Installing from GitHub will allow you to install the latest, albeit unstable, version of `dysh`.
To install the main branch of `dysh` with all extra dependencies from GitHub:

.. code::

    $ pip install "dysh[all] @ git+https://github.com/GreenBankObservatory/dysh"
