********************************
Testing & GitHub Integrations
********************************

Test Code
=========

All new features **must** have accompanying test code that verifies their correctness; no PR will be accepted without it.  A jupyter notebook explaining how to use the new feature is also strongly recommended.   Bug fixes should have a regression test added to verify the issue is fixed and doesn't return.  Tests must be runnable by pytest.   If a test writes a file to disk, it **must** use`pytest's `tmp_path  fixture <https://docs.pytest.org/en/latest/how-to/tmp_path.html>`_ to ensure the file is removed after the test.

Tests must be named `test_*.py` and located within a `tests` directory of their parent module. All tests are run as a GitHub Workflow in dysh's continuous integration suite defined in `dysh/.github/workflows/ci.yml`.


Pre-Commit Hooks
================

We have several pre-commit hooks that ensure committed code follows desired standards. They can be found in `dysh/.pre-commit-config.yaml`.   pre-commit is automatically installed if you used the `dev` dependency group to install dysh (`uv sync --dev`).  `pre-commit` will run whenever you commit code.

.. code-block:: bash

    (dysh) $ git add newfile.py
    (dysh) $ git commit -m "adding new python file" newfile.py

The `pre-commit` will run  to make sure all of the staged files are formatted correctly. You'll see a message like this, then the information about your commit.

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
