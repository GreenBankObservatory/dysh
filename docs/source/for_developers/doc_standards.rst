*************
Documentation
*************

Creating Documentation
======================

All Python classes and methods should be documented and the documentation must follow `the NumPy convention. <https://numpydoc.readthedocs.io/en/latest/format.html>`_.  Use of major features should be demonstrated in an example notebook.

Mermaid Diagrams
================

Diagrams can be directly added to the documentation by using the `sphinxcontrib-mermaid <https://github.com/mgaitan/sphinxcontrib-mermaid>`_ package.
Here's an example:

.. mermaid::

    flowchart LR
        A[Item 1] --> B[Item 2]
        B --> C[Item 3]

To learn more, see the `package documentation <https://sphinxcontrib-mermaid-demo.readthedocs.io/en/latest/>`_.
Mermaid also offers an `online editor <https://mermaid.live>`_ which can be used to design diagrams.

Building documentation
======================

.. code::

    $ cd docs/source
    $ uv run make html

Sphinx Autobuilds
=================

Sphinx autobuilds allow live checking of local documentation.
The docs will be published at ``http://127.0.0.1:8000/``;  this can be modified with arguments ``--host <hostname> --port <port number>``.

.. code::

    $ uv run sphinx-autobuild docs/source src/dysh/docs/build -b html

    The HTML pages are in docs/build.
    [sphinx-autobuild] Serving on http://127.0.0.1:8000
    [sphinx-autobuild] Waiting to detect changes...


Any changes in the ``dysh/docs/`` directory can be viewed at the output URL, and will auto-update when changes are made.
Changes to the API (i.e. ``src/dysh``) won't be automatically updated.
To close the server, simply ``CTRL+C``.
