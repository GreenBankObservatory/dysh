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
from messages import TerminalMessages
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
from dysh.fits.sdfitsload import SDFITSLoad #, get_size, Obsblock, baseline
from dysh.spectra.obsblock import Obsblock

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

        self.main_text = QLabel("Welcome to the Dysh GUI")
        self.button = QPushButton("Select file")
        self.button.clicked.connect(self.get_files)
        self.main_layout.addWidget(self.main_text, 0, 0, 1, 1)
        self.main_layout.addWidget(self.button, 1, 0, 1, 1)
        
        self.setCentralWidget(self.main_widget)

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
    def __init__(self, fpath):
        """ Initializes the main window """
        super(MainWindow, self).__init__()
        self.fpath = fpath
        #self.fpath = "/Users/victoriacatlett/Desktop/TGBT21A_501_11.raw.vegas.fits"

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

    def _load_data(self):
        """ Opens up the FITS file """
        # [TODO] Load lists in a QThread so the main screen can be created
        self.fdata = SDFITSLoad(self.fpath)
        self.fdata._loadlists(1)

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
    terminal_messages = TerminalMessages()
    terminal_messages.welcome()
    
    app = QApplication(sys.argv)
    apply_stylesheet(app, theme='dark_purple.xml')
    win = StartWindow()

    sys.exit(app.exec_())
    
    # [TODO] Add functions to execute on exit
    # terminal_messages.goodbye()

if __name__ == '__main__':
    start()
