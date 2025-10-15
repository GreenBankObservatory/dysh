*******************
GitHub Integrations
*******************

Actions
=======

Tests should be named `test_*.py` and located within a `tests` directory of their parent module. All tests are run in dysh's continuous integration suite defined in `dysh/.github/workflows/ci.yml`.


Pre-Commit Hooks
================

We have several pre-commit hooks that ensure committed code follows desired standards. They can be found in `dysh/.pre-commit-config.yaml`.   pre-commit is automattically installed if you used the `dev` dependency group to install dysh (`uv sync --dev`).  `pre-commit` will run whenever you commit code.

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

If any check Failed, then a file was modified by pre-commit, so you must commit again to get those changes.

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
