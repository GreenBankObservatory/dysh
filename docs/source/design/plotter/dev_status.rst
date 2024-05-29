*********************
Status of Development
*********************

I've been having some trouble with using just `matplotlib` for the interactive plotter. 

PyQt with PyQtGraph
===================

One option is to add the `PyQt` dependency back in. 

Pros:
  * `PyQtGraph` is based off of `matplotlib`
  * Allows future extensibility into a full GUI

Cons:
  * Large package download


Plotly with PyQt
================

`Plotly` makes a lot of this interactive stuff easy. 

Pros:
  * Lots of built-in interactivity options
  * Easy to export interactive plots into Sphinx

Cons:
  * Not based off of `matplotlib`
  * Default behavior is to create an HTML file which gets rendered by a browser. Can specify renderer.
  * Using `PyQt`` to create a popup window introduces the `PyQt`` dependency