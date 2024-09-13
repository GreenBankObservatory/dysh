*****************
Design and Status
*****************

The plotter is based off of `PyQt` with `PyQtGraph`.

Plot Types
==========
Each type of plot will be defined by its own plot class with a specific `PyQtGraph` widget.

`SingleSpectrum`
--------------

A `SingleSpectrum` plot will display a single frquency-vs-intensity `spectrum` object.

Waterfall
---------

This plot type has not been implemented.

Selection Areas
===============

Different plot types may require different types of selection.

Single-Point selection
----------------------

On spectra, users may one to be able to click on individual points to create a selection. This feature has not been implemented yet.

Rectangular Region of Interest (ROI)
------------------------------------

A selection object on the plot is defined as a rectangular  Region of Interest (ROI).
