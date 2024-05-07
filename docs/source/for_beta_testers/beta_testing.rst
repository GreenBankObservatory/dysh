*****************************
Instructions for Beta Testers
*****************************

For beta testing we ask that your provide feedback in the form of:

* Reporting issues. For example, if something does not work or the results are not accurate.
* Making suggestions. For example, requesting new features.
* Letting us know what you think.

If you are interested in beta testing `dysh`, please contact Pedro Salas (psalas@nrao.edu) to be added to the beta testers list. Thanks for your interest!


Providing feedback
==================

If you encounter a problem with `dysh`, would like to request a new feature or enhancement or would like to leave feedback, please do so using `GitHub issues <https://github.com/GreenBankObservatory/dysh/issues>`_. There are some basic instructions of how to do this :ref:`here <githubissues>`. This requires `creating a free account <https://github.com/>`_ on GitHub if you do not have one.

If you prefer not to create a GitHub account, please provide your feedback to the `dysh-beta mailing list <https://groups.google.com/g/dysh-beta/about>`_, or send an email to dysh-beta@googlegroups.com. Additionally, we will provide a `form for collecting feedback <https://forms.gle/27tg9adfLbDnUyz37>`_.

When providing feedback, please provide

* `Python` version
* `dysh` version
* operating system

If reporting an issue please also provide

* the input data (either as a link or a location inside the GBO computing environment)
* a minimum working example to reproduce the error
* outputs (for example, error messages or figures)
* any additional information that might help us reproduce and understand the problem


Example feedback
----------------

Here are examples of feedback on GitHub

* reporting an issue, `Issue #88 <https://github.com/GreenBankObservatory/dysh/issues/88>`_
* requesting a modification, `Issue #78 <https://github.com/GreenBankObservatory/dysh/issues/78>`_

.. _beta-install:

Installing `dysh`
=================

Here we provide additional installation steps that include creating a virtual environment to keep `dysh` isolated from your system `Python` version.
We provide steps for working in one of `GBO data reduction hosts <https://greenbankobservatory.org/science/gbt-observers/public-access-data-reduction/>`_ (e.g., fourier), and if you're working outside one of the GBO data reduction hosts.

.. tab:: At GBO

    Create a `Python3.11` virtual environment

    .. code:: bash

        /users/gbosdd/python/bin/python3.11 -m venv /home/scratch/$USER/dysh-0.2-env

    Activate your virtual environment

    .. code:: bash

        source /home/scratch/$USER/dysh-0.2-env/bin/activate

    Install `dysh` and Jupyter Lab

    .. code:: bash

        pip install jupyterlab dysh==0.2.0b

.. tab:: Outside of GBO

    Create a `Python3.9+` virtual environment

    .. code:: bash

        python3 -m venv /path/to/venv

    Activate your virtual environment

    .. code:: bash

        source /path/to/venv/bin/activate

    Install `dysh` and Jupyter Lab

    .. code:: bash

        pip install jupyterlab dysh==0.2.0b

In the future we will provide `dysh` executables at GBO.
