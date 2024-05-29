*******************
Choosing a Renderer
*******************

There are several different situations in which one might use the plotter, and it must render and be interactive in all of them. 

Scripting
=========

If a user wants the interactive plotter to appear during or after a script execution, it must appear as a pop-up window. 

iPython
=======

The plotter must appear when called during an iPython session, and it must stay open and respond while new commands are entered. 

Notebooks
=========

The plotter must render beneath a notebook cell when it is called. The plot need only update when the cell is re-run. 