*************
Documentation
*************

Setting Up Sphinx Autobuilds
============================

Here are the steps to set up Sphinx autobuilds so that you can check your documentation edits live.

1. First, navigate to your `dysh` root directory and activate the `hatch` environment.

.. code-block:: bash

    $ hatch shell

2. Next, tell hatch to run the docs. The docs will be published at `http://127.0.0.1:8000/`.

.. code-block:: bash

    (dysh) $ hatch run docs

3. If you would like the docs to publish at a specific host and port, such as `http://thales:9876`, then add the appropriate flags:

.. code-block:: bash

    (dysh) $ hatch run docs --host thales --port 9876

4. You may now make changes in the `dysh/docs/` directory and see the live changes at the appropriate URL in your browser. To close the server, simply `CTRL+C`.

Docstring Format
================

All Python functions must contain a docstring which follows the NumPy convention. You can learn more about this convention here: https://numpydoc.readthedocs.io/en/latest/format.html

Mermaid Diagrams
================

Diagrams can be directly in these text files by using the `sphinxcontrib-mermaid` package. Here's an example:

.. mermaid::

    flowchart LR
        A[Item 1] --> B[Item 2]
        B --> C[Item 3]

To learn more, see the `package documentation <https://sphinxcontrib-mermaid-demo.readthedocs.io/en/latest/>`_. Mermaid also offers an `online editor <https://mermaid.live>`_ which can be used to design diagrams.
