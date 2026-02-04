"""
Script for preparing notebooks so they render nicely in Read the Docs.
"""

import nbformat as nbf


def inject_static_gui(filename):
    """ """

    # Load the notebook.
    nb = nbf.read(filename, nbf.NO_CONVERT)

    # The code to inject.
    cmd = """from dysh.plot.staticgui import StaticGUI
from dysh import plot
plot.switch_frontend(StaticGUI)
%matplotlib inline"""

    # Create a new cell and add a tag so it doesn't show.
    new_cell = nbf.v4.new_code_cell(cmd)
    new_cell["metadata"]["tags"] = "remove-cell"

    # Add the cell at the start of the notebook.
    nb["cells"].insert(0, new_cell)

    # Write the notebook.
    with open(filename, "w") as f:
        nbf.write(nb, f)
