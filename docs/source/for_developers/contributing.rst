*************************
How to Contribute to Dysh
*************************

tl;dr: Create an issue if you find a bug. Create a pull request if you want to contribute code

Setting up a Dysh Development Environment
=========================================

In order to contribute to Dysh, you will need the following:

In order to set up a Dysh development environment, you will need the following:
#. `git <https://git-scm.com/>`_
#. `Python 3.9 <https://www.python.org/downloads/>`_
#. `Hatch <https://hatch.pypa.io/latest/install/>`_
#. `Pre-commit <https://pre-commit.com/#install>`_


See {TODO} for help on setting up these pre-requisites on your development platform
Then, you can:

.. code-block:: bash

    git clone https://github.com/GreenBankObservatory/dysh.git
    cd dysh
    pre-commit install  # Install the pre-commit hooks into your git repo
    hatch shell  # Create virtual env, install dependencies, and activate the environment


Contributing Code
=================

All code contributions will require:
#. `A GitHub account <https://github.com/signup>`_
#. `A fork the Dysh repo <https://github.com/GreenBankObservatory/dysh/fork>`_
#. An issue to be created


Two example workflows:

Feature Request/Bug Report
++++++++++++++++++++++++++

1. Check the `Dysh Issues <https://github.com/GreenBankObservatory/dysh/issues?q=>`_ to see if the bug/feature has already been reported/requested
2. If not, `create an Issue <https://github.com/GreenBankObservatory/dysh/issues/new>`_
3. Dysh developers will communicate with you via the Issue to help diagnose the bug

Bugs will be fixed as developer time allows.

I want to add a feature/fix a bug
+++++++++++++++++++++++++++++++++

I want to fix a bug
If you think you know how to fix the bug, we would welcome a bug fix!

1. In (a clone of) your fork, create a "feature branch" with a descriptive name like ``desc-of-your-bug``
2. Commit code to fix the bug
3. Create a `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_
4. Wait for Dysh developers to communicate with you further!


6. Check the `Dysh Issues <https://github.com/GreenBankObservatory/dysh/issues?q=>`_ to see if a similar feature has already been requested
7. If not, `create an Issue <https://github.com/GreenBankObservatory/dysh/issues/new>`_ describing the feature
8. Dysh developers will communicate with you via the Issue to help diagnose the bug
