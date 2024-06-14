# PACKAGE IMPORTS
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from pyqtgraph import PlotWidget
from pyqtgraph.jupyter import GraphicsLayoutWidget
import pyqtgraph as pg
import sys #, os, psutil, getpass, socket
# import wget
#import numpy as np
#import pyqtgraph as pg
#from astropy.io import fits
#from time import time
#import pandas as pd
#import argparse
from screeninfo import get_monitors
from qt_material import apply_stylesheet
from IPython.display import display

# LOCAL GUI IMPORTS
# from widgets.tables import FITSHeaderTable
# from widgets.graphs import SingleSpectrum
# from widgets.QIPython import QIPythonConsoleWidget
# from util.dataload import FITSFileDialog

# from util.core import ThreadCallbacks, DyshWorker
# from widgets.splash import SplashScreen
# from widgets.graphs import *
# from widgets.layouts import *
# from util.dataload import DataLoader

# DYSH IMPORTS
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.plot.renderer import Renderer

# PARALLELIZATION
# from concurrent.futures import ThreadPoolExecutor
# from threading import Thread

class CustomGLW(GraphicsLayoutWidget):
    def get_frame(self):
        # rather than eating up cpu cycles by perpetually updating "Updating plot",
        # we will only update it opportunistically on a redraw.
        # self.request_draw()
        #update()
        return super().get_frame()
    
class SingleSpectrum(PlotWidget):
    """ Spectrum Plot """
    def __init__(self,spectrum,**kwargs):
        super().__init__()
        self._spectrum = spectrum
        self.add_data()

    def config(self):
        # [TODO] connect spec_num to the hline value
        spec_num = 0
        self.setLabels(left='Intensity', bottom='Frequency')
        self.setTitle(f"Spectrum {spec_num}")

    def update_data(self, spectrum):
        self.clear()
        self._spectrum = spectrum
        self.add_data()

    def add_data(self, **kwargs):
        #self._plot_kwargs.update(kwargs)
        #self.get_kwargs(**kwargs)

        s = self._spectrum
        sa = s.spectral_axis
        sf = s.flux

        self.plot(sa, sf)

    def get_kwargs(self, **kwargs):
        self.lw =  self._plot_kwargs['linewidth']
        self.xunit = self._plot_kwargs["xaxis_unit"]
        self.yunit = self._plot_kwargs["yaxis_unit"]

class iSingleSpectrum:
    """ Spectrum Plot """
    def __init__(self,spectrum,**kwargs):
        self.app = pg.mkQApp()
        self.win = CustomGLW(css_width="1000px", css_height="600px")
        self._spectrum = spectrum
        self.add_data()

    def config(self):
        # [TODO] connect spec_num to the hline value
        spec_num = 0
        self.setLabels(left='Intensity', bottom='Frequency')
        self.setTitle(f"Spectrum {spec_num}")

    def update_data(self, spectrum):
        self.clear()
        self._spectrum = spectrum
        self.add_data()

    def add_data(self, **kwargs):
        #self._plot_kwargs.update(kwargs)
        #self.get_kwargs(**kwargs)

        s = self._spectrum
        sa = s.spectral_axis
        sf = s.flux

        self.myplot = self.win.addPlot(title="My plot")
        plot_indx = self._spectrum.mask == False
        self.myplot.plot(sa[plot_indx], sf[plot_indx], color="white")
        self.myplot.plot(sa[~plot_indx], sf[~plot_indx], color="red")

    def get_kwargs(self, **kwargs):
        self.lw =  self._plot_kwargs['linewidth']
        self.xunit = self._plot_kwargs["xaxis_unit"]
        self.yunit = self._plot_kwargs["yaxis_unit"]

class DyshMainWindow(QMainWindow):
    """The main window of the GUI"""

    def __init__(self, spectrum, fpath=None):
        """Initializes the main window"""
        super(DyshMainWindow, self).__init__()

        self.spectrum = spectrum

        self.setWindowTitle("Dysh GUI")
        self._init_geometry(0.8)

        self.info_threads()
        self._init_main_panel()
        
        self.show()

    def _init_geometry(self, mult):
        """
        Draws the GUI on the primary monitor

        Parameters
        ----------
            mult : int or float
                proportion of total size to draw window (0.8 = 80%)

        """
        for m in get_monitors():
            if m.is_primary:
                self.width = int(m.width * mult)
                self.height = int(m.height * mult)
                self.xpos = int(m.x + (m.width * (1 - mult)) / 2)
                self.ypos = int(m.y + (m.height * (1 - mult)) / 2)
        self.setGeometry(self.xpos, self.ypos, self.width, self.height)

    def info_threads(self):
        """Updates info on available threads"""
        self.threadCountActive = QThreadPool.globalInstance().activeThreadCount()
        self.threadCountTotal = QThreadPool.globalInstance().maxThreadCount()
        # print(f"You are using {self.threadCountActive} of the {self.threadCountTotal} available QThreads")

    def _init_main_panel(self):
        #self._clear_all()
        #self._load_data()
        self._init_UI()

    #@SingleThread
    def SDFITS_load_all(self, fpath):
        self.sdfits = GBTFITSLoad(fpath)

    def _load_data(self):
        """Opens up the FITS file"""
        # [TODO] Load lists in a QThread so the main screen can be created
        # [TODO] Add logic to determine if GBTFITSLoad or another
        #s_load = DyshWorker(target=self.SDFITS_load_all, args=(self.fpath, 1))
        #s_load.start()
        #url = "https://www.gb.nrao.edu/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
        #self.fpath = wget.download(url)
        self.fpath = "TGBT21A_501_11.raw.vegas.fits"

        self.SDFITS_load_all(self.fpath) #s_load.join()
        self.scan = self.sdfits.getps(152, ifnum=0, plnum=0)
        self.scan.calibrate()
        self.fdata = self.scan.timeaverage(weights='tsys')

    def _init_UI(self):
        """Creates the skeleton structure of the GUI"""
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)
        self._init_plots()
        self._init_ROI()
        self.main_layout.addWidget(self.spec_plot, 0, 0, 1, 1)
        self._init_select_button()

    def _init_select_button(self):
        self.select_button = QPushButton("Save Selection")
        self.main_layout.addWidget(self.select_button, 1, 0, 1, 1)
        self.select_button.clicked.connect(self.get_ROI)

    def _init_ROI(self):
        self.roi = pg.RectROI((1.39E9, -1), (1E6, 2), pen=pg.mkPen('r', width=2), movable=True, resizable=True)
        self.spec_plot.addItem(self.roi)

    def get_ROI(self):
        #self.roi.

    def _init_plots(self):
        """Creates the plot canvases"""
        # [TODO] Do this in a QThread so the main screen can be created
        self.spec_plot = SingleSpectrum(self.spectrum)
        # print(f"NINTEGRATIONS: {self.fdata.nintegrations(1)}")
        # self.waterfall = WaterfallSpectrum(self.fdata)
        #self.tab3_layout.addWidget(self.spec_plot, 0, 0, 1, 2)

    def _init_plot_sidebar(self):
        pass

    def _clear_all(self):
        while self.main_layout.count():
            child = self.main_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def closeEvent(self, *args, **kwargs):
        super(QMainWindow, self).closeEvent(*args, **kwargs)
        #self.terminal.stop()

class App(QApplication):
    def __init__(self, spectrum, *args):
        QApplication.__init__(self, *args)
        self.main = DyshMainWindow(spectrum)
        self.main.show()

def makePlot(spectrum):
    rend = Renderer()
    print("RENDERER: ", rend.render_type)
    if rend.render_type == "notebook":
        plot = iSingleSpectrum(spectrum)
        plot.show()
    else:
        app = App(spectrum, sys.argv)
        apply_stylesheet(app, theme="dark_purple.xml")
        app.exec_()

if __name__ == "__main__":
    from dysh.fits.gbtfitsload import GBTFITSLoad

    # filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
    filename = "E:/Code/GitHub/Work/forks/dysh/test_file.fits"
    sdfits = GBTFITSLoad(filename)
    psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
    ta = psscan.timeaverage(weights="tsys")
    # ta.mask[0:3000] = True
    # breakpoint()

    makePlot(ta)