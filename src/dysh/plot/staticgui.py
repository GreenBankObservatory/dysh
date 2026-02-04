"""
dysh headless "GUI".
"""

from matplotlib.backends.backend_agg import FigureCanvasAgg

from .basegui import BaseGUI


class StaticGUI(BaseGUI):
    """
    Agg based GUI, i.e., no GUI at all.
    """

    def __init__(self, plotbase):
        self.canvas = FigureCanvasAgg(plotbase.figure)
