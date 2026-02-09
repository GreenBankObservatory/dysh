"""
dysh headless "GUI".
"""

from matplotlib.backends.backend_agg import FigureCanvasAgg


class StaticGUI:
    def __init__(self, plotbase):
        self.canvas = FigureCanvasAgg(plotbase.figure)

    def show(self):
        return
