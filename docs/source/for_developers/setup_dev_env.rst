************************************
Setting Up a Development Environment
************************************


Pre-Requisites
==============

There are a few core things you will need to develop Dysh:

1. A Linux computer (Mac and Windows should also work, but we don't have docs for those)
2. Python (specifically you should be using the minimum-supported version of Python; see ``pyproject.toml``)
3. `pre-commit <https://pre-commit.com/>`_: used to enforce code quality checks
4. `Hatch <https://hatch.pypa.io/latest/install/#pipx>`_ (optional, technically)

Python
------

See site-specific instructions [TODO] for GBO, UMD, etc. instructions -- the necessary version of Python is probably already available somewhere.


If you do not already have access to the minimum-supported version of Python, you will need to install it. This can be done manually via `the official site <https://www.python.org/downloads/>`_, but it is much easier to use pyenv.


Once you can run ``python --version`` and see the correct version pop up, you are ready to move on.

pre-commit
----------



Dysh Development Environment
============================

Start to finish, how do you initialize a Dysh development environment?


Via Hatch
---------

.. code-block:: sh

    git clone git@github.com:GreenBankObservatory/dysh.git
    cd dysh
    hatch shell  #
    hatch run test  # Prove that the tests can run and pass

Manual Virtual Environment
--------------------------

In some environments it makes more sense to manage your environments yourself. For example, at GBO we have a variety of hosts that all mount a common filesystem, but can have different OS versions. Here the standard "one env per project per user" breaks down, and it is better to manage things yourself.

.. code-block:: sh

    git clone git@github.com:GreenBankObservatory/dysh.git
    cd dysh
    python3.9 -m venv /path/to/venv  # Create your venv
    source /path/to/venv/bin/activate  # Activate your venv
    pip install -e .[all]  # Install dysh in edit mode, including all extra deps
    hatch run test  # Use hatch to trigger pytest inside your manually-created venv


Advanced Cases
==============

If you need to exactly replicate the Python dependencies of either:
1. A given CI deployment
2. A given Dysh deployment at GBO (or elsewhere)

You will instead be relying on the lockfile (``requirements.txt``) to provide exact version pins of all of your dependencies.

With the same version of Python as the environment you are trying to replicate, and with the same OS (or as similar as possible):

.. code-block:: sh

    python3.9 -m venv /path/to/venv  # Create your venv
    source /path/to/venv/bin/activate  # Activate your venv
    pip install -r requirements.txt  # Install from lockfile
