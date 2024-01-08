import getpass
import os
import socket
import sys

import psutil
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from pyqtgraph import GraphicsLayoutWidget, ImageItem


class FITSHeaderTable(QWidget):
    """Table of FITS Header information"""

    def __init__(self):
        """Initializes the table widget"""
        super().__init__()
        self.make_layout()

    def make_layout(self):
        self.title = QLabel("FITS Header")
        self.tbl = QTableWidget()
        self.tbl_layout = QVBoxLayout()
        self.setLayout(self.tbl_layout)
        self.tbl_layout.addWidget(self.title)
        self.tbl_layout.addWidget(self.tbl)

    def load(self, data):
        """
        Gets the keys

        Parameters
        ----------
            data : dict
                A dictionary of the FITS header

        """
        ks = [k for k in data.keys()]

        self.tbl.setRowCount(len(ks))
        self.tbl.setColumnCount(2)
        self.tbl.setHorizontalHeaderLabels(["Header Key", "Header Value"])

        for i, ki in enumerate(ks):
            self.tbl.setItem(i, 0, QTableWidgetItem(str(ki)))
            self.tbl.setItem(i, 1, QTableWidgetItem(str(data[ki])))


class FITSDataTable(QTableWidget):
    """Table of FITS Header information"""

    def __init__(self):
        """Initializes the table widget"""
        super().__init__()

    def make_layout(self):
        self.title = QLabel("FITS Data")
        self.tbl = QTableWidget()
        self.tbl_layout = QVBoxLayout()
        self.setLayout(self.tbl_layout)
        self.tbl_layout.addWidget(self.title)
        self.tbl_layout.addWidget(self.tbl)

    def get_keys(self, data):
        """
        Gets the keys

        Parameters
        ----------
            data : dict
                A dictionary of the FITS column names

        """
        ks = data.keys()

        self.tbl.setRowCount(len(ks))
        self.tbl.setColumnCount(4)
        self.tbl.setHorizontalHeaderLabels(["Header Key", "Value", "Unit", "TFORM"])

        for i, ki in enumerate(ks):
            try:
                self.tbl.setItem(i, 0, QTableWidgetItem(str(ki)))
                self.tbl.setItem(i, 1, QTableWidgetItem(str(self.hdr_df[ki][0])))
                self.tbl.setItem(i, 2, QTableWidgetItem(str(self.h_data_info[ki]["TUNIT"])))
                self.tbl.setItem(i, 3, QTableWidgetItem(str(self.h_data_info[ki]["TFORM"])))
            except:
                print(f"Issue encountered for {ki} (data)")
