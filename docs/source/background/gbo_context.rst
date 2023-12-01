**************************
Context Within GBO Systems
**************************

Background
==========

GBTIDL
------

`GBTIDL` is the current software used to calibrate SDFITS files from the GBT.

Requirements
============

At its core, `dysh` needs to read in an SDFITS file and ultimately output calibrated data.

.. mermaid::

    flowchart LR
        A[SDFITS File] --> B[Dysh]
        B --> C[Output Data]
