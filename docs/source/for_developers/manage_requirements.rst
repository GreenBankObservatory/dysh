*********************
Managing Requirements
*********************

Add any new package dependency to `pyproject.toml`. If the package is needed for `dysh` functionality, add it to the `dependencies = [...]` list.
If the package is only needed for development, add it to the `dev = [...]` list. If the package is needed for notebooks, add it to the `nb = [...]` list.  The current lists in `pyproject.toml` are listed below.   Do not remove any of these without prior discussion with the `dysh` maintainers.

**Do not edit  requirements.txt.** It is managed by `uv`.

Once the changed files files  are to GitHub, verify that all GitHub Actions tests pass to ensure that these requirements work for all supported versions of Python and all supported operating systems.

.. code-block:: Python

 # excerpted from pyproject.toml

    dependencies = [
      "httpx",
      "astropy>=6.1",
      "astroquery",
      "ipython",
      "jplephem",
      "matplotlib",
      "numpy<2",
      "pandas",
      "rich",
      "scipy",
      "specutils>=2",
      "tenacity",
    ]

    [project.optional-dependencies]
    nb = [
      "jupyter",
      "jupyterlab",
      "ipympl",
      "jupyter-nbextensions-configurator>=0.6.4",
    ]

    all = ["dysh[nb]"]

    [dependency-groups]
    dev = [
      "coverage[toml]",
      "ipdb",
      "myst-nb",
      "nbclient",
      "nbformat",
      "numpydoc",
      "pre-commit",
      "pytest",
      "pytest-cov",
      "pytest-xdist>=3.6.1",
      "ruff",
      "sphinx",
      "sphinx-autobuild",
      "sphinx-book-theme",
      "sphinx-copybutton",
      "sphinx-design",
      "sphinx-inline-tabs",
      "sphinxcontrib-mermaid",
      "dysh[all]", # ensure that all extras are included in dev envs
      "pytest-timeout>=2.4.0",
      "yamlfmt>=1.1.1",
   ]
