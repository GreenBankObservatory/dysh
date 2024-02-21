*******************
Installing ``dysh``
*******************

``dysh`` requires Python 3.9+ and recent versions of
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

.. warning::
    `dysh` is currently in development and the above command will install the latest stable version of `dysh` which might not reflect the contents of the documentation. For beta testing please see :ref:`beta-install`.

From github
===========

To install from github without creating a separate virtual environment:

.. code::

    $ git clone git@github.com:GreenBankObservatory/dysh.git
    $ cd dysh
    $ pip install -e .

If you wish to install using a virtual environment, which we strongly recommend if you plan to contribute to the code, see :doc:`installation instructions for developers <../for_developers/install>`.
