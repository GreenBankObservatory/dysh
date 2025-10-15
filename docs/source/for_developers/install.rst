************************************
Installation and Virtual Environment
************************************

For developing `dysh` code,
we recommend the use of a python virtual environment for development. The example below uses `uv <https://docs.astral.sh/uv/>`_.

Using `uv`
----------

#. Clone the repo, install `uv`, and create the virtual environment, by default in dysh/.venv:

    .. code-block:: bash

        $ git clone git@github.com:GreenBankObservatory/dysh.git
        $ pip install uv
        $ cd dysh
        $ uv venv

#. Install dysh and its required packages. This will create a uv.lock file containing the exact dependencies.

    .. code-block:: bash

       $ uv sync --all-extras --dev

#. Verify the installation

    .. code-block:: bash

       # should point to $CWD/.venv/bin/dysh
       $ uv run which dysh

       # should print dysh version and exit
       $ uv run dysh --version

#. Optional: run the test suite

    .. code-block:: bash

       $ uv run pytest -n auto

#. Optional: build the documentation locally

    .. code:: bash

      $ cd docs/source
      $ uv run make html
      $ xdg-open _build/html/index.html
