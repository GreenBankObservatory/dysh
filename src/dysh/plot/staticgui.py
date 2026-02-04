"""
dysh headless "GUI".
"""

from IPython.display import display
from matplotlib.backends.backend_agg import FigureCanvasAgg

from .basegui import BaseGUI


class StaticGUI(BaseGUI):
    def __init__(self, plotbase):
        self.canvas = FigureCanvasAgg(plotbase.figure)

    def show(self):
        display(self.canvas.figure)
