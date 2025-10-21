*********************
Managing Requirements
*********************

New package dependencies  must be added with `uv add`.  If the package is needed for `dysh` functionality:

.. code:: bash

 $ uv add <package name>

If the package is only needed for development, add it to the `dev` group:

.. code:: bash

 $ uv add --dev <package name>

If the package is needed for notebooks, add it to the `nb` optional dependencies group.

.. code:: bash

   $ uv add --group nb <package name>

`uv add` will update uv.lock and pyproject.toml; you should not edit those files directly.

Do not remove any packages without prior discussion with the `dysh` maintainers.

Once the changed files files are pushed to GitHub, verify that all GitHub Actions tests pass to ensure that these requirements work for all supported versions of Python and all supported operating systems.
