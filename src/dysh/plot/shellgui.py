"""
Graphical User Interface for dysh shell.
"""

import tkinter as tk
from tkinter import ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from ..util.core import in_notebook

if not in_notebook():
    root = tk.Tk()
    root.withdraw()


class ShellGUI:
    def __init__(self, plotbase):
        self.root = tk.Toplevel()
        self.root.geometry("800x800")
        self.root.title("dysh")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        self.frame = ttk.Frame(self.root, width=800, height=500)
        self.frame.pack(expand=True, fill=tk.BOTH, side=tk.TOP)
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=1)
        self.frame.rowconfigure(0, weight=1)
        self.frame.rowconfigure(1, weight=1)
        # This combination seems to work well through FastX.
        plotbase.figure.set_dpi(100)
        plotbase.figure.set_figwidth(9)
        plotbase.figure.set_figheight(6)
        self.canvas = FigureCanvasTkAgg(plotbase.figure, master=self.frame)
        # Using a manager causes things to freeze!
        # self.manager = FigureManagerTk(self.canvas, 0, self.root)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root, pack_toolbar=False)
        self.toolbar.update()
        self.canvas.draw()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X, expand=True)
        self.canvas.get_tk_widget().grid(column=0, row=1, columnspan=2, sticky="nsew")

        if plotbase._selector is not None:
            self.button_clear = ttk.Button(
                master=self.frame, text="Clear All Regions", command=plotbase._selector.clear_regions
            )
            self.button_clear.grid(column=0, row=0)
            self.button_clear_one = ttk.Button(
                master=self.frame, text="Clear Region", command=plotbase._selector.clear_region
            )
            self.button_clear_one.grid(column=1, row=0)

    def show(self):
        return
