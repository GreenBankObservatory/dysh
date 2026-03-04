"""
dysh lab (Jupyter) plotting interface.
"""

import functools
import warnings

from ipympl.backend_nbagg import Canvas, FigureManager
from IPython.display import display
from ipywidgets import Button, Dropdown, HBox, VBox
from matplotlib.backends.backend_agg import FigureCanvasAgg

from .basegui import BaseGUI


class LabGUI(BaseGUI):
    """
    NBAgg based GUI.
    Used for jupyter.
    """

    def __init__(self, plotbase):
        # Suppress traitlets DeprecationWarning from ipympl Toolbar MRO issue (ipympl#488).
        # warnings.filterwarnings alone is unreliable in Jupyter kernels, so use
        # catch_warnings for a guaranteed local suppression.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning, module="traitlets")
            self.canvas = Canvas(plotbase.figure)
            self.canvas.manager = FigureManager(self.canvas, 0)
        plotbase.figure.canvas.header_visible = False  # Remove canvas header.

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
        self.buttonbar = [self.button_clear, self.button_clear_one]

        if hasattr(plotbase, "set_xaxis_unit"):
            self.dropdown_xunit_opts = [(unit, i) for i, unit in enumerate(plotbase._xaxis_unit_options())]
            self.dropdown_xunit = Dropdown(
                options=self.dropdown_xunit_opts, value=0, description="x-axis units:", tooltip="Set x-axis units"
            )
            self.buttonbar.append(self.dropdown_xunit)
        if hasattr(plotbase, "set_yaxis_unit"):
            self.dropdown_yunit_opts = [(unit, i) for i, unit in enumerate(plotbase._yaxis_unit_options())]
            self.dropdown_yunit = Dropdown(
                options=self.dropdown_yunit_opts, value=0, description="y-axis units:", tooltip="Set y-axis units"
            )
            self.buttonbar.append(self.dropdown_yunit)

        self.ui = VBox([HBox(self.buttonbar), self.canvas])

    def connect_buttons(self, plotbase):
        if plotbase.has_selector():
            # Button widget takes care of not registering the same callback twice.
            self.button_clear.on_click(plotbase._selector.clear_regions)
            self.button_clear_one.on_click(plotbase._selector.clear_region)
        if hasattr(plotbase, "set_xaxis_unit"):
            self.dropdown_xunit.observe(
                functools.partial(self._dropdown_xunit, dropdown_func=plotbase.set_xaxis_unit),
                names="value",
                type="change",
            )
        if hasattr(plotbase, "set_yaxis_unit"):
            self.dropdown_yunit.observe(
                functools.partial(self._dropdown_yunit, dropdown_func=plotbase.set_yaxis_unit),
                names="value",
                type="change",
            )

    def disconnect_buttons(self):
        self.button_clear.on_click(self.button_clear._click_handlers.callbacks[0], remove=True)
        self.button_clear_one.on_click(self.button_clear_one._click_handlers.callbacks[0], remove=True)

    def is_window_alive(self):
        return True

    def show(self):
        display(self.ui)

    def _dropdown_xunit(self, change, dropdown_func=None):
        opt_dict = {v: k for k, v in self.dropdown_xunit_opts}
        dropdown_func(opt_dict[change["new"]])

    def _dropdown_yunit(self, change, dropdown_func=None):
        opt_dict = {v: k for k, v in self.dropdown_yunit_opts}
        dropdown_func(opt_dict[change["new"]])


class StaticLabGUI(BaseGUI):
    """
    Agg based GUI, i.e., no GUI at all.
    Used for jupyter rendering in documentation.
    """

    def __init__(self, plotbase):
        self.canvas = FigureCanvasAgg(plotbase.figure)

    def show(self):
        display(self.canvas.figure)
