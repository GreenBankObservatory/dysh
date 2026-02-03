"""
dysh lab (Jupyter) plotting interface.
"""

from ipympl.backend_nbagg import Canvas, FigureManager
from IPython.display import display
from ipywidgets import Button, HBox, VBox


class LabGUI:
    def __init__(self, plotbase):
        # figure.set_dpi(100)
        plotbase.figure.set_figwidth(10)
        plotbase.figure.set_figheight(6)
        self.canvas = Canvas(plotbase.figure)
        self.canvas.manager = FigureManager(self.canvas, 0)
        plotbase.figure.canvas.header_visible = False  # Remove canvas header.

        if hasattr(plotbase, "_selector"):
            self.button_clear = Button(
                description="Clear All Regions",
                disabled=False,
                tooltip="Clear all regions from plot",
            )
            self.button_clear.on_click(plotbase._selector.clear_regions)
            self.button_clear_one = Button(
                description="Clear Region",
                disabled=False,
                tooltip="Clear selected region from plot",
            )
            self.button_clear_one.on_click(plotbase._selector.clear_region)

            self.ui = VBox([HBox([self.button_clear, self.button_clear_one]), plotbase.figure.canvas])

        else:
            self.ui = VBox([plotbase.figure.canvas])

    def show(self):
        display(self.ui)
