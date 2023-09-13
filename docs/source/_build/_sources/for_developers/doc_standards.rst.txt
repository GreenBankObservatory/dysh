***********************
Documentation Standards
***********************

Sphinx Autobuilds
=====================

* First, you need to create a documentation environment.

.. code-block:: bash

    cd dysh/docs
    /path/to/python -m venv docs-env
    source docs-env/bin/activate
    pip install -r requirements.txt

* Next, set up the build command. 

.. code-block:: bash

    cp .env.template .env

* Add values for `DOCS_ROOT`, `DOCS_HOST`, and `DOCS_PORT` in `.env`
* Start the autobuild

.. code-block:: bash
    
    source .env
    startdocs

* Go to `http://{$DOCS_HOST}:{$DOCS_PORT}` in a web browser. You should now see the documentation with live edits as you save changes. 

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