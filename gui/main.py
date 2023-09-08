# PACKAGE IMPORTS
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import sys, os, psutil, getpass, socket
import numpy as np
import pyqtgraph as pg
from astropy.io import fits
from time import time
import pandas as pd
import argparse
from screeninfo import get_monitors
from qt_material import apply_stylesheet

# ADD PATHS
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('.'))

# LOCAL GUI IMPORTS
from widgets.tables import FITSHeaderTable
from widgets.graphs import *
from dataload import FITSFileDialog
'''
from widgets.splash import SplashScreen
from widgets.graphs import *
from gui.core import Worker
from gui.dataload import DataLoader
'''
# DYSH IMPORTS
from dysh.util.messages import *
from dysh.fits.sdfitsload_parallel import SDFITSLoad #, get_size, Obsblock, baseline
from dysh.fits.gbtfitsload import GBTFITSLoad #, get_size, Obsblock, baseline
from dysh.spectra.obsblock import Obsblock

# PARALLELIZATION
from concurrent.futures import ThreadPoolExecutor
from threading import Thread

class ThreadCallbacks:
    def progress(future):
        pass
        #print('.', end='', flush=True)

class DyshWorker(QObject, Thread):
    """ Thread to run a function that returns a value """
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}):
        self.finished = pyqtSignal()
        self.progress = pyqtSignal(int)
        Thread.__init__(self, group=group, target=target, name=name, args=args, kwargs=kwargs)
        self._return = None
        #print(self.__dir__())
    
    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)
    
    def join(self):
        Thread.join(self)
        return self._return
    
class StartWindow(QMainWindow):
    """ The startup window of the GUI """
    def __init__(self):
        """ Initializes the startup window """
        super().__init__()
        self.setWindowTitle("Dysh GUI")
        self._init_geometry(0.8)
        self._init_UI()
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
                self.xpos = int(m.x + (m.width * (1-mult))/2)
                self.ypos = int(m.y + (m.height * (1-mult))/2)
        self.setGeometry(self.xpos, self.ypos, self.width, self.height)

    def _init_UI(self):
        """ Creates the skeleton structure of the GUI """
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

        # Make the UI Items
        self.main_text = QLabel("Welcome to the Dysh GUI")
        self.button = QPushButton("Select file")
        self.button.clicked.connect(self.get_files)

        self._init_site_selection()

        # Add the UI Items
        self.main_layout.addWidget(self.main_text, 0, 0, 1, 2)
        self.main_layout.addWidget(self.combo_telescope, 1, 0, 1, 1)
        self.main_layout.addWidget(self.combo_rx, 1, 1, 1, 1)
        self.main_layout.addWidget(self.button, 2, 0, 1, 2)
        
        self.setCentralWidget(self.main_widget)

    def _init_site_selection(self):
        self.combo_telescope = QComboBox()
        self.combo_telescope.addItem('Auto-Detect')
        self.combo_telescope.addItem('Green Bank Telescope (GBT)')
        self.combo_telescope.addItem('Green Bank 20-meter Telescope')
        self.combo_telescope.addItem('Large Millimeter Telescope (LMT)')

        self.combo_rx = QComboBox()
        self.update_combo_rx()

        self.combo_telescope.currentIndexChanged.connect(self.update_combo_rx)

    def update_combo_rx(self):
        # [TODO] Load the RX info from a JSON file
        self.combo_rx.clear()
        self.combo_rx.setEnabled(True)

        ci = int(self.combo_telescope.currentIndex())
        if ci == 0:
            # Auto-Detect
            self.combo_rx.setEnabled(False)
        elif ci == 1:
            # Green Bank Telescope (GBT)
            self.combo_rx.addItem("PF1 (0.29 - 0.395 GHz)")
            self.combo_rx.addItem("PF1 (0.385 - 0.52 GHz)")
            self.combo_rx.addItem("PF1 (0.51 - 0.69 GHz)")
            self.combo_rx.addItem("PF1 (0.68 - 0.92 GHz)")
            self.combo_rx.addItem("PF2 (0.9 - 1.23 GHz)")
            self.combo_rx.addItem("L (1.15 - 1.73 GHz)")
            self.combo_rx.addItem("S (1.73 - 2.6 GHz)")
            self.combo_rx.addItem("UWBR (0.7 - 4.0 GHz)")
            self.combo_rx.addItem("C (3.95 - 8.0 GHz)")
            self.combo_rx.addItem("X (8.0 - 12.0 GHz)")
            self.combo_rx.addItem("Ku (12.0 - 15.4 GHz)")
            self.combo_rx.addItem("KFPA (17.0 - 27.5 GHz)")
            self.combo_rx.addItem("Ka F1 (26.0 - 31.0 GHz)")
            self.combo_rx.addItem("Ka F2 (30.5 - 37.0 GHz)")
            self.combo_rx.addItem("Ka F3 (36.0 - 39.5 GHz)")
            self.combo_rx.addItem("Q (38.2 - 49.8 GHz)")
            self.combo_rx.addItem("W1 (68.0 - 74.0 GHz)")
            self.combo_rx.addItem("W2 (73.0 - 80.0 GHz)")
            self.combo_rx.addItem("W3 (79.0 - 86.0 GHz)")
            self.combo_rx.addItem("W4 (85.0 - 92.0 GHz)")
            self.combo_rx.addItem("ARGUS (75.0 - 115.5 GHz)")
        elif ci == 2:
            # Green Bank 20-meter telescope
            self.combo_rx.addItem("L (1.15 - 1.73 GHz)")
            self.combo_rx.addItem("X (8.0 - 12.0 GHz)")
        elif ci == 3:
            # Green Bank 20-meter telescope
            self.combo_rx.addItem("RSR")
            self.combo_rx.addItem("SEQUOIA")
            self.combo_rx.addItem("MSIP1mm")
            self.combo_rx.addItem("B4R")
            self.combo_rx.addItem("TolTEC")

    def get_files(self):
        # [TODO] Figure out why this makes you do it twice?
        file_dialog = FITSFileDialog()
        if file_dialog.exec_() == QDialog.Accepted:
            self.fpath = file_dialog.selectedFiles()[0]
            self.open_main_window()

    def open_main_window(self):
        self.close()
        self.new_window = MainWindow(self.fpath)
        self.new_window.show()
        

class MainWindow(QMainWindow):
    """ The main window of the GUI """
    def __init__(self, fpath=None):
        """ Initializes the main window """
        super(MainWindow, self).__init__()
        self.fpath = "/Users/victoriacatlett/Desktop/TGBT21A_501_11.raw.vegas.fits"#fpath

        self.setWindowTitle("Dysh GUI")
        self._init_geometry(0.8)

        self.info_threads()
        self._load_data()
        self._init_UI()
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
                self.xpos = int(m.x + (m.width * (1-mult))/2)
                self.ypos = int(m.y + (m.height * (1-mult))/2)
        self.setGeometry(self.xpos, self.ypos, self.width, self.height)

    def info_threads(self):
        """ Updates info on available threads """
        self.threadCountActive = QThreadPool.globalInstance().activeThreadCount()
        self.threadCountTotal = QThreadPool.globalInstance().maxThreadCount()
        print(f"You are using {self.threadCountActive} of the {self.threadCountTotal} available QThreads")

    def SDFITS_load_all(self, fpath, n):
        fdata = SDFITSLoad(fpath)
        fdata._loadlists(n)
        return fdata

    def _load_data(self):
        """ Opens up the FITS file """
        # [TODO] Load lists in a QThread so the main screen can be created
        # [TODO] Add logic to determine if GBTFITSLoad or another
        s_load = DyshWorker(target=self.SDFITS_load_all, args=(self.fpath,1))
        s_load.start()
        self.f_data = s_load.join()

    def _init_UI(self):
        """ Creates the skeleton structure of the GUI """
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

        self._init_tables()
        self._init_plots()

    def _init_tables(self):
        """ Creates tables of FITS information """
        # [TODO] Add selection logic for if len(bintable) > 1
        # [TODO] Do this in a QThread so the main screen can be created
        self.hdr0_tbl = FITSHeaderTable()
        self.hdr0_tbl.load(self.fdata.primaryheader())
        self.hdr1_tbl = FITSHeaderTable()
        self.hdr1_tbl.load(self.fdata.binheader()[0])
        self.main_layout.addWidget(self.hdr0_tbl, 0, 0, 1, 1)
        self.main_layout.addWidget(self.hdr1_tbl, 0, 1, 1, 1)

    def _init_plots(self):
        """ Creates the plot canvases """
        # [TODO] Do this in a QThread so the main screen can be created
        self.fdata.summary()
        #print(f"NINTEGRATIONS: {self.fdata.nintegrations(1)}")
        #self.waterfall = WaterfallSpectrum(self.fdata)



def start():
    friendly_messages = FriendlyMessages()
    friendly_messages.welcome()
    
    app = QApplication(sys.argv)
    apply_stylesheet(app, theme='dark_purple.xml')
    #win = StartWindow()
    win = MainWindow()
    sys.exit(app.exec_())
    
    # [TODO] Add functions to execute on exit
    # friendly_messages.goodbye()

if __name__ == '__main__':
    start()
