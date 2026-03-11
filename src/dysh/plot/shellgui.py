"""
Graphical User Interface for dysh shell.
"""

import functools
import os

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from ..util.core import in_notebook
from .basegui import BaseGUI

no_display = os.environ.get("DISPLAY", "") == ""
if not no_display:
    import tkinter as tk
    from tkinter import ttk

if not in_notebook() and not no_display:
    root = tk.Tk()
    root.withdraw()


class ShellGUI(BaseGUI):
    """
    Tk based GUI.
    Used for plotting in the shell.
    """

    def __init__(self, plotbase):
        self.root = tk.Toplevel()
        self.root.geometry("900x700")
        self.root.title("dysh")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        self.bframe = ttk.Frame(self.root, width=900, height=100)  # Button frame.
        self.bframe.pack(expand=True, fill=tk.BOTH, side=tk.TOP)
        self.frame = ttk.Frame(self.root, width=900, height=300)
        self.frame.pack(expand=True, fill=tk.BOTH, side=tk.TOP)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)
        self.canvas = FigureCanvasTkAgg(plotbase.figure, master=self.frame)
        # Using a manager causes things to freeze!
        # self.manager = FigureManagerTk(self.canvas, 0, self.root)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root, pack_toolbar=False)
        self.toolbar.update()
        self.canvas.draw()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X, expand=True)
        self.canvas.get_tk_widget().grid(column=0, row=0, sticky="nsew")

        # Buttons.
        self.button_clear = ttk.Button(master=self.bframe, text="Clear All Regions")
        self.button_clear.grid(column=0, row=0, sticky="NS")
        self.button_clear_one = ttk.Button(master=self.bframe, text="Clear Region")
        self.button_clear_one.grid(column=1, row=0)

        # Dropdown menus.
        if hasattr(plotbase, "set_xaxis_unit"):
            self.xunit = tk.StringVar(self.root)
            options = plotbase._xaxis_unit_options()
            self.dropdown_xunit = ttk.OptionMenu(
                self.bframe,
                self.xunit,
                options[0],
                *options,
                command=functools.partial(self._dropdown_xunit, plotbase=plotbase),
            )
            self.dropdown_xunit.config(width=8)
            self.dropdown_xunit.grid(column=2, row=0, sticky="NS")
        if hasattr(plotbase, "set_yaxis_unit"):
            self.yunit = tk.StringVar(self.root)
            options = plotbase._yaxis_unit_options()
            self.dropdown_yunit = ttk.OptionMenu(
                self.bframe,
                self.yunit,
                options[0],
                *options,
                command=functools.partial(self._dropdown_yunit, plotbase=plotbase),
            )
            self.dropdown_yunit.grid(column=3, row=0, sticky="NS")

    def connect_buttons(self, plotbase):
        if plotbase.has_selector():
            self.button_clear.bind("<Button>", plotbase._selector.clear_regions)
            self.button_clear_one.bind("<Button>", plotbase._selector.clear_region)

    def disconnect_buttons(self):
        self.button_clear.unbind("<Button>")
        self.button_clear_one.unbind("<Button>")

    def is_window_alive(self):
        return bool(self.root.winfo_exists())

    def show(self):
        return

    def _dropdown_xunit(self, event=None, plotbase=None):
        if hasattr(plotbase, "set_xaxis_unit"):
            plotbase.set_xaxis_unit(self.xunit.get())

    def _dropdown_yunit(self, event=None, plotbase=None):
        if hasattr(plotbase, "set_yaxis_unit"):
            plotbase.set_yaxis_unit(self.yunit.get())
