# PACKAGE IMPORTS
import sys  # , os, psutil, getpass, socket

# PARALLELIZATION
from concurrent.futures import ThreadPoolExecutor
from threading import Thread

import wget
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from qt_material import apply_stylesheet

# import numpy as np
# import pyqtgraph as pg
# from astropy.io import fits
# from time import time
# import pandas as pd
# import argparse
from screeninfo import get_monitors
from util.core import DyshWorker, ThreadCallbacks
from util.dataload import DataLoader, FITSFileDialog
from widgets.graphs import *
from widgets.layouts import *
from widgets.QIPython import QIPythonConsoleWidget
from widgets.splash import SplashScreen

# LOCAL GUI IMPORTS
from widgets.tables import FITSHeaderTable

from dysh.fits.gbtfitsload import GBTFITSLoad

# DYSH IMPORTS
from dysh.util.messages import *
from dysh.util.parallelization import SingleThread


class SelectPanel(QGridLayout):
    """The startup window of the GUI"""

    def __init__(self):
        """Initializes the startup window"""
        super().__init__()
        self._init_UI()

    def _init_UI(self):
        """Creates the skeleton structure of the GUI"""

        # Make the UI Items
        self.main_text = QLabel("Welcome to the Dysh GUI")
        self.button = QPushButton("Select file")
        self.button.clicked.connect(self.get_files)

        self._init_site_selection()

        # Add the UI Items
        self.addWidget(self.main_text, 0, 0, 1, 2)
        self.addWidget(self.combo_telescope, 1, 0, 1, 1)
        self.addWidget(self.combo_rx, 1, 1, 1, 1)
        self.addWidget(self.button, 2, 0, 1, 2)

    def _init_site_selection(self):
        self.combo_telescope = QComboBox()
        self.combo_telescope.addItem("Auto-Detect")
        self.combo_telescope.addItem("Green Bank Telescope (GBT)")
        self.combo_telescope.addItem("Green Bank 20-meter Telescope")
        self.combo_telescope.addItem("Large Millimeter Telescope (LMT)")

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
        self.file_dialog = FITSFileDialog()


class DyshMainWindow(QMainWindow):
    """The main window of the GUI"""

    def __init__(self, fpath=None):
        """Initializes the main window"""
        super(DyshMainWindow, self).__init__()
        FriendlyMessages.hello()

        self.setWindowTitle("Dysh GUI")
        self._init_geometry(0.8)

        self.info_threads()
        # self._init_select_panel()
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

    def _init_select_panel(self):
        self.main_widget = QWidget()
        self.main_layout = SelectPanel()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)
        if self.main_layout.file_dialog.exec_() == QDialog.Accepted:
            self.fpath = self.main_layout.file_dialog.selectedFiles()[0]

    def _init_main_panel(self):
        # self._clear_all()
        self._load_data()
        self._init_UI()

    # @SingleThread
    def SDFITS_load_all(self, fpath):
        self.sdfits = GBTFITSLoad(fpath)

    def _load_data(self):
        """Opens up the FITS file"""
        # [TODO] Load lists in a QThread so the main screen can be created
        # [TODO] Add logic to determine if GBTFITSLoad or another
        # s_load = DyshWorker(target=self.SDFITS_load_all, args=(self.fpath, 1))
        # s_load.start()
        # url = "https://www.gb.nrao.edu/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
        # self.fpath = wget.download(url)
        self.fpath = "TGBT21A_501_11.raw.vegas.fits"

        self.SDFITS_load_all(self.fpath)  # s_load.join()
        self.scan = self.sdfits.getps(152, ifnum=0, plnum=0)
        self.scan.calibrate()
        self.fdata = self.scan.timeaverage(weights="tsys")

    def _init_UI(self):
        """Creates the skeleton structure of the GUI"""
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

        self._init_sidebar()
        self._init_toggle_btn()
        self._init_tabs()

        self.main_layout.addWidget(self.toggle_btn, 0, 0, 1, 1)
        self.main_layout.addWidget(self.sidebar, 1, 0, 1, 1)

        self.main_layout.addWidget(self.tabs, 0, 1, 2, 2)
        self._init_tables()
        self._init_plots()
        self._init_terminal()

    def _init_tabs(self):
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()

        self.tab1_layout = QGridLayout()
        self.tab2_layout = QGridLayout()
        self.tab3_layout = QGridLayout()
        self.tab4_layout = QGridLayout()

        self.tab1.setLayout(self.tab1_layout)
        self.tab2.setLayout(self.tab2_layout)
        self.tab3.setLayout(self.tab3_layout)
        self.tab4.setLayout(self.tab4_layout)

        self.tabs.addTab(self.tab1, "File")
        self.tabs.addTab(self.tab2, "Waterfall")
        self.tabs.addTab(self.tab3, "Calibrated Spectrum")
        self.tabs.addTab(self.tab4, "Console")

    def _init_sidebar(self):
        self.sidebar = CollapsibleSideBar()
        self.sidebar.add_box(title="my box 1", contentWidget=QLabel("content 1"))
        self.sidebar.add_box(title="my box 2", contentWidget=QLabel("content 2"))

    def _init_toggle_btn(self):
        self.toggle_btn = QPushButton("Dock")
        self.toggle_btn.clicked.connect(self.toggle_hidden)

    def toggle_hidden(self):
        if self.sidebar.isHidden() == True:
            self.sidebar.setHidden(False)
        else:
            self.sidebar.setHidden(True)

    def _init_tables(self):
        """Creates tables of FITS information"""
        # [TODO] Add selection logic for if len(bintable) > 1
        # [TODO] Do this in a QThread so the main screen can be created
        self.hdr0_tbl = FITSHeaderTable()
        self.hdr0_tbl.load(self.sdfits.primaryheader())
        self.hdr1_tbl = FITSHeaderTable()
        self.hdr1_tbl.load(self.sdfits.binheader()[0])
        self.tab1_layout.addWidget(self.hdr0_tbl, 0, 0, 1, 1)
        self.tab1_layout.addWidget(self.hdr1_tbl, 0, 1, 1, 1)

    def _init_plots(self):
        """Creates the plot canvases"""
        # [TODO] Do this in a QThread so the main screen can be created
        self.spec_plot = SingleSpectrum(self.fdata)
        # print(f"NINTEGRATIONS: {self.fdata.nintegrations(1)}")
        # self.waterfall = WaterfallSpectrum(self.fdata)
        self.tab3_layout.addWidget(self.spec_plot, 0, 0, 1, 2)

    def _init_plot_sidebar(self):
        pass

    def _init_terminal(self):
        self.terminal = QIPythonConsoleWidget()
        self.tab4_layout.addWidget(self.terminal, 0, 0, 1, 1)

    def _clear_all(self):
        while self.main_layout.count():
            child = self.main_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def closeEvent(self, *args, **kwargs):
        super(QMainWindow, self).closeEvent(*args, **kwargs)
        self.terminal.stop()
        FriendlyMessages.goodbye()


class App(QApplication):
    def __init__(self, *args):
        QApplication.__init__(self, *args)
        self.main = DyshMainWindow()
        self.main.show()


def main(args):
    # global app
    app = App(args)
    apply_stylesheet(app, theme="dark_purple.xml")
    app.exec_()


if __name__ == "__main__":
    main(sys.argv)
