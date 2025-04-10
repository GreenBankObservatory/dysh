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

If you prefer not to create a GitHub account, please provide your feedback to the `dysh-beta mailing list <https://groups.google.com/g/dysh-beta/about>`_, or send an email to dysh-beta@googlegroups.com. Additionally, we provide a `form for collecting feedback <https://forms.gle/gf9rydgNE8v7iDKR8>`_.

When providing feedback, please indicate

* `Python` version (:code:`python --version`)
* `dysh` version (:code:`dysh --version`)
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
We provide steps for working in one of `GBO data reduction hosts <ihttps://greenbankobservatory.org/portal/gbt/processing/#data-reduction-machines>`_ (e.g., fourier), and if you're working outside one of the GBO data reduction hosts. The provided command will install `dysh` and `jupyterlab`, which can be used to run the example notebooks locally.
Please remember that as of version 0.4 `dysh` requires `Python3.10+`.

.. tab:: At GBO

    On a linux data reduction machine using a terminal.

    Create a `Python3.11` virtual environment

    .. code:: bash

        /users/gbosdd/python/bin/python3.11 -m venv /home/scratch/$USER/dysh-env

    Activate your virtual environment

    .. code:: bash

        source /home/scratch/$USER/dysh-env/bin/activate

    Install `dysh` and Jupyter Lab

    .. code:: bash

        pip install "dysh[nb]==0.4.2"

    Check that the correct version of `dysh` was installed

    .. code:: bash

        dysh --version

    It should print 0.4.2.

.. tab:: Outside of GBO

    From a terminal.

    Create a `Python3.10+` virtual environment

    .. code:: bash

        python3 -m venv /path/to/venv

    Activate your virtual environment

    .. code:: bash

        source /path/to/venv/bin/activate

    Install `dysh` and Jupyter Lab

    .. code:: bash

        pip install "dysh[nb]==0.4.2"

    Check that the correct version of `dysh` was installed

    .. code:: bash

        dysh --version

    It should print 0.4.2.

`dysh` is installed in the GBO data reduction hosts, however, it may not be the latest version.
You can launch it using

.. code:: bash

    dysh


Previous beta releases
======================

Feedback on previous beta releases is also welcome. Here you can find links to previous beta release documents.

0.3.0
-----

`Instructions <https://docs.google.com/document/d/182FMM3f0pi54r6qDc_Ttv59Sgms2i4hJzVqf9Nw9GfY/edit?usp=sharing>`_ and `questionnaire <https://forms.gle/MGSD2tR1sPdZXxNq9>`_.

0.2.0
-----

`Instructions <https://docs.google.com/document/d/1RrHaiwmrDnPbMLdNY99_hBZzyWyYKsw0UCM8FKqhIKo/edit?usp=sharing>`_ and `questionnaire <https://forms.gle/27tg9adfLbDnUyz37>`_.
