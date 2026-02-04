"""
dysh lab (Jupyter) plotting interface.
"""

from ipympl.backend_nbagg import Canvas, FigureManager
from IPython.display import display
from ipywidgets import Button, HBox, VBox

from .basegui import BaseGUI


class LabGUI(BaseGUI):
    def __init__(self, plotbase):
        self.canvas = Canvas(plotbase.figure)
        self.canvas.manager = FigureManager(self.canvas, 0)
        self.canvas.header_visible = False  # Remove canvas header.

        self.button_clear = Button(
            description="Clear All Regions",
            disabled=False,
            tooltip="Clear all regions from plot",
        )
        self.button_clear_one = Button(
            description="Clear Region",
            disabled=False,
            tooltip="Clear selected region from plot",
        )

        self.ui = VBox([HBox([self.button_clear, self.button_clear_one]), self.canvas])

    def connect_buttons(self, plotbase):
        if plotbase.has_selector():
            # Button widget takes care of not registering the same callback twice.
            self.button_clear.on_click(plotbase._selector.clear_regions)
            self.button_clear_one.on_click(plotbase._selector.clear_region)

    def disconnect_buttons(self):
        self.button_clear.on_click(self.button_clear._click_handlers.callbacks[0], remove=True)
        self.button_clear_one.on_click(self.button_clear_one._click_handlers.callbacks[0], remove=True)

    def is_window_alive(self):
        return True

    def show(self):
        display(self.ui)
