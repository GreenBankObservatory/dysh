*************
Documentation
*************

Setting Up Sphinx Autobuilds
============================

Here are the steps to set up Sphinx autobuilds so that you can check your documentation edits live. 

1. First, navigate to your `dysh` root directory and activate the `hatch` environment.

.. code-block:: bash

    $ hatch shell

2. Next, copy the environment file template.

.. code-block:: bash

    (dysh) $ cp .env.template .env

3. Add values for ``DOCS_ROOT``, ``DOCS_HOST``, and ``DOCS_PORT`` in `.env`
4. Start the autobuild

.. code-block:: bash
    
    (dysh) $ source .env
    (dysh) $ startdocs

5. Go to `http://{$DOCS_HOST}:{$DOCS_PORT}` in a web browser. You should now see the documentation with live edits as you save changes. 

.. note::
    Do not commit the `.env` file to `git`. 

Docstring Format
================

Gotta format the docstrings

Mermaid Diagrams
================

Diagrams can be directly in these text files by using the `sphinxcontrib-mermaid` package. Here's an example:

.. mermaid::

    flowchart LR
        A[Item 1] --> B[Item 2]
        B --> C[Item 3]

To learn more, see the `package documentation <https://sphinxcontrib-mermaid-demo.readthedocs.io/en/latest/>`_. Mermaid also offers an `online editor <https://mermaid.live>`_ which can be used to design diagrams. 