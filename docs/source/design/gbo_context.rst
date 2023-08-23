**************************
Context Within GBO Systems
**************************

Overview 
========

At its core, `dysh` needs to read in an SDFITS file and ultimately output calibrated data.

.. mermaid::

    flowchart LR
        A[SDFITS File] --> B[Dysh]
        B --> C[Output Data]