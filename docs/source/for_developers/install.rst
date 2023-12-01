***************************
Installation for Developers
***************************

Here are the steps to install ``dysh`` if you want to develop code.

We use `hatch <https://hatch.pypa.io/>`_ to manage the build environment.
The usual caveats apply how you set up your python development environment.

#. Clone the repo and install hatch.

    .. code-block:: bash

        $ git clone git@github.com:GreenBankObservatory/dysh.git
        $ cd dysh
        $ pip install hatch

#. Create and activate a virtual environment with hatch and install the packages required for development.
   The virtual environment will be created the first time; subsequent invoking ``hatch shell`` will simply load the created environment.

    .. code-block:: bash

        $ hatch shell
        (dysh) $ pip install -r requirements_dev.txt

#. Build and install the package

    .. code-block:: bash

        (dysh) $ hatch build
        (dysh) $ pip install -e .

#. You can exit this environment (which effectively had started a new shell).

    .. code-block:: bash

        (dysh) $ exit
        $

#. Each time when you come back in this directory without being in this virtual environment, you'll need to load the virtual environment.

    .. code-block:: bash

        $ hatch shell

Notice you can *only* do that from within the original install directory tree.

Additional Installation Options
-------------------------------

.. include:: install_developer.rst

Testing
=======
We use `pytest` for unit and integration testing.
From the top-level dysh directory, run:

.. code-block:: bash

    $ pytest
