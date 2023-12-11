*************
Documentation
*************

How does it work, and how can you contribute to it?

How to Contribute to Dysh Docs
==============================

Local Builds
++++++++++++

The first thing you'll want to do is get the docs building locally. This makes it very easy to develop, and you can be *nearly certain* that the docs will build on Read the Docs if they have successfully built locally (assuming you've rememembered to commit all your files...)

To run a live docs server (i.e. reload on changes)

.. code-block:: sh

    cd "$DYSH_REPO_ROOT"
    hatch run docs

You can pass any arguments through to ``sphinx-autobuild``, e.g. to bind to a specific host and port:

.. code-block:: sh

    cd "$DYSH_REPO_ROOT"
    hatch run docs --host 0.0.0.0 --port 8001



Flowcharts/Diagrams
+++++++++++++++++++

Where possible, please prefer Mermaid to pre-rendered flowcharts/diagrams. It makes it easy to embed the full source in the docs, which makes future edits significantly easier.

Dysh supports Mermaid Diagrams via the ``sphinxcontrib-mermaid`` package. Here's an example:

.. mermaid::

    flowchart LR
        A[Item 1] --> B[Item 2]
        B --> C[Item 3]

To learn more, see the `sphinxcontrib-mermaid documentation <https://sphinxcontrib-mermaid-demo.readthedocs.io/en/latest/>`_. Mermaid also offers an `online editor <https://mermaid.live>`_ which can be used to design diagrams.

Read the Docs
=============

Dysh's documentation is hosted on Read the Docs, at https://readthedocs.org/projects/dysh/. Docs are built directly from the Dysh GitHub repo.

RTD hosts several different active versions of Dysh documentation:

- ``latest``: built from the latest commit on ``main``
- ``stable``: built from that latest stable version tag
- ``release-X.Y``/``vX.Y``: built from the latest commit on a given release branch
