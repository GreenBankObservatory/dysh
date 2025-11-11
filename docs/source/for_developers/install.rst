************************************
Installation and Virtual Environment
************************************

For developing `dysh` code,
we recommend the use of a python virtual environment. The example below uses `uv <https://docs.astral.sh/uv/>`_.   Before installing `dysh`, developers should install `uv` following one of the methods in the `uv docs <https://docs.astral.sh/uv/getting-started/installation/>`_.  `uv` is the only tool that can sync the environment to the lockfile, so to install the known working development environment, `uv` is needed.

Using `uv`
----------

#. Clone the repo, install ``dysh`` and its required packages. This will use the uv.lock file containing the exact dependencies as well as a virtual environment directory in `dysh/.venv`.

    .. code-block:: bash

       $ git clone git@github.com:GreenBankObservatory/dysh.git
       $ cd dysh
       $ uv sync --all-extras

#. Verify the installation.  Below, `uv run` is used to execute commands, however one could also run `source .venv/bin/activate` which eliminates the need to prepend `uv run` to each command.

    .. code-block:: bash

       # should point to $CWD/.venv/bin/dysh
       $ uv run which dysh

       # should print dysh version and exit
       $ uv run dysh --version

#. Optional: run the test suite

    .. code-block:: bash

       $ uv run pytest

  You can optionally add `-n auto` to spawn multiple runners for pytest.

#. Optional: build the documentation locally

    .. code:: bash

      $ cd docs/source
      $ uv run make html
      $ xdg-open _build/html/index.html
