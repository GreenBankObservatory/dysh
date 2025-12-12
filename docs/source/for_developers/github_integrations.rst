********************************
Testing & GitHub Integrations
********************************

Test Code
=========

All new features **must** have accompanying test code that verifies their correctness; no PR will be accepted without it.  A jupyter notebook explaining how to use the new feature is also strongly recommended.   Bug fixes should have a regression test added to verify the issue is fixed and doesn't return.  Tests must be runnable by pytest.   If a test writes a file to disk, it **must** use`pytest's `tmp_path  fixture <https://docs.pytest.org/en/latest/how-to/tmp_path.html>`_ to ensure the file is removed after the test.   The `dysh` development team can provide help developing regression tests for bug fixes.

Tests must be named `test_*.py` and located within a `tests` directory of their parent module. All tests are run as a GitHub Workflow in dysh's continuous integration suite defined in `dysh/.github/workflows/ci.yml`.


Pre-Commit Hooks
================

We use `pre-commit <https://pre-commit.com/>` with
hooks that ensure committed code follows desired standards. The hooks be found in `dysh/.pre-commit-config.yaml`.

There two steps for enabling pre-commit in your repo: 1) installing the pre-commit package and 2) installing the hooks in your clone of the dysh repo.

Installing the pre-commit package
-----------------------------------------------

You can install `pre-commit` globally, so that it is active for any of your projects that have hooks, or just locally for `dysh`.  For global install, in your normal environment:

.. code-block:: bash

    $ uv tool install pre-commit

If you prefer local to `dysh`,  the `pre-commit` package is automatically installed with `dev` dependency group when installing `dysh`:

.. code-block:: bash

    (dysh) $ uv sync --dev

Installing dysh's pre-commit hooks
------------------------------------------------
To install the pre-commit hooks into your repo (`.git/hooks/pre-commit`)

Globally:

.. code-block:: bash

   pre-commit install

Locally:

.. code-block:: bash

   (dysh) $uv run pre-commit install


`pre-commit` will run  whenever you commit to make sure all of the staged files are formatted correctly.

.. code-block:: bash

    (dysh) $ git add newfile.py
    (dysh) $ git commit -m "adding new python file" newfile.py

You'll see a message like this, then the information about your commit.

.. code-block:: bash
   :emphasize-lines: 11

    trim trailing whitespace.................................................Passed
    fix end of files.........................................................Passed
    check toml...........................................(no files to check)Skipped
    check for added large files..............................................Passed
    debug statements (python)................................................Passed
    detect private key.......................................................Passed
    mixed line ending........................................................Passed
    check docstring is first.................................................Passed
    check for case conflicts.................................................Passed
    yamlfmt..............................................(no files to check)Skipped
    ruff format.............................................................. Failed

    - hook id: ruff-format
    - files were modified by this hook

    1 file reformatted

     ruff check...............................................................Passed

If a check Failed, then a file was modified by pre-commit, so you must commit again to get those changes.

.. code-block:: bash

    $  git commit -m "adding new python file" newfile.py

    trim trailing whitespace.................................................Passed
    fix end of files.........................................................Passed
    check toml...........................................(no files to check)Skipped
    check for added large files..............................................Passed
    debug statements (python)................................................Passed
    detect private key.......................................................Passed
    mixed line ending........................................................Passed
    check docstring is first.................................................Passed
    check for case conflicts.................................................Passed
    yamlfmt..............................................(no files to check)Skipped
    ruff format..............................................................Passed
    ruff check...............................................................Passed
    [your_branch_name  commit_hash] "adding new python file"
     1 file changed, 2 insertions(+), 2 deletions(-)
