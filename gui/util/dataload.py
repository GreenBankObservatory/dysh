import argparse
import os
import sys
from time import time

import numpy as np
import pandas as pd
from PyQt5.QtCore import *

# from PyQt5.QtCore import QRunnable, QObject, Qt, QThreadPool, pyqtSignal, pyqtSlot, QMutex
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# sys.path.insert(0, os.path.abspath('../src'))
from dysh.fits.sdfitsload import SDFITSLoad


class FITSFileDialog(QFileDialog):
    """File dialog to select a FITS file"""

    def __init__(self):
        """Initializes the file dialog widget"""
        super().__init__()
        self.get_fpaths()
        self.config()
        self.getOpenFileName(caption="Open an SDFITS file")

    def config(self):
        self.setDirectory(self.f_root)
        # self.setFilter(self.f_filter)

    def get_fpaths(self):
        # [TODO] Have this check what system you're on and build the QFileSystemModel
        # [TODO] Add the FITS file filter
        self.f_root = QDir("/home/dysh/example_data/")
        # self.f_filter = "FITS Files (*.fits)"
        # self.f_model = QFileSystemModel()
        # self.f_model.setRootPath("/home/dysh/example_data/")

    def get_selection(self):
        # [TODO] Closing the file dialog without selecting anything produces an error. Fix that
        fdialog = QFileDialog(self, "Select an SDFITS File", froot, filter_sdfits)
        if fdialog.exec_() == QDialog.Accepted:
            fpath = fdialog.selectedFiles()[0]


class DataLoader(QObject):
    finished = pyqtSignal()
    my_data = pyqtSignal()
    my_hdr = ()

    def load_data(self, fpath):
        print("RUNNABLE: LOADING DATA")
        t0 = time()
        p0, a0 = self.get_memory_usage()

        parser = argparse.ArgumentParser(prog="revisedstructure")
        parser.add_argument("--wcs", action="store_true", help="create WCS for each Spectrum1D")
        parser.add_argument("--fix", action="store_true", help="fix warnings about keywords in WCS creation")
        parser.add_argument("--profile", "-p", action="store_true", help="run the profiler")
        parser.add_argument("--baseline", "-b", action="store_true", help="remove baselines")
        parser.add_argument(
            "--maxload", "-m", action="store", help="maximum number of spectra to load (create obsblocks)", default=1e16
        )
        args = parser.parse_args()

        data_sl = SDFITSLoad(fpath)
        for h in range(1, len(data_sl._hdu)):
            data_sl._loadlists(fix=args.fix, wcs=args.wcs, hdu=int(h), maxspect=float(args.maxload))

        t1 = time()
        p1, a1 = self.get_memory_usage()

        for o in data_sl._obsblock:
            data = []
            for oi in range(len(o)):
                data.append(o[oi].spectral_axis.data)
            # self.data = np.asarray(data).T
            # print(f'NEW DATA: {np.shape(self.data)}, {type(self.data)}')
            print(f"SHAPE: {np.shape(o)}, {np.shape(o[0].spectral_axis)}")
            print(o)
        print(f"It took {t1-t0:.2f} seconds to load the file with the SDFITS Loader.")
        print(f"This process added {p1-p0:.2f} MiB of RAM usage.\n")

    def load_header(self, fdata):
        print("RUNNABLE: LOADING HEADER")
        # [TODO] These are only loading one value, but they may change for each scan
        all_hks = []
        self.fhdr = {}
        self.tfhdr = {}
        # self.fhdr_d = {}
        for hk in self.fdata[1].header:
            hkv = str(fd[1].header[hk])
            all_hks.append(hkv)

        n_tfield = int(fdata[1].header["TFIELDS"])
        n_points = int(fdata[1].header["NAXIS1"])

        # PLACES TO SAVE THE INFO
        self.h_fits = {}
        self.h_data_info = {}
        h_data_all = {}

        # GET INFO ON THE WHOLE FILE
        h_fits_info = [
            "XTENSION",
            "BITPIX",
            "NAXIS",
            "NAXIS1",
            "NAXIS2",
            "PCOUNT",
            "GCOUNT",
            "TFIELDS",
            "TELESCOP",
            "BACKEND",
            "SITELONG",
            "SITELAT",
            "SITEELEV",
            "EXTNAME",
        ]
        for hfi in h_fits_info:
            try:
                self.h_fits[hfi] = fd[1].header[hfi]
                # print(f'Header key {hfi} has value {fd[1].header[hfi]}')
                # all_hks.remove(hfi)
            except:
                print(f"Issue with header {hfi}")

        # GET SCAN-SPECIFIC INFORMATION
        for i in range(n_tfield):
            # Use i+1 because the TFIELDs use 1-based counting, while range() uses 0-based
            try:
                tf0 = fdata[1].header[f"TTYPE{i+1}"]
                tf_data = fdata[1].data[tf0]
                all_hks.remove(tf0)
            except:
                tf0 = np.zeros(n_points)
            try:
                tf1 = fdata[1].header[f"TFORM{i+1}"]
                all_hks.remove(tf1)
            except:
                tf1 = "-"
            try:
                tf2 = fdata[1].header[f"TUNIT{i+1}"]
                all_hks.remove(tf2)
            except:
                tf2 = "-"
            self.h_data_info[tf0] = {"TFORM": tf1, "TUNIT": tf2}
            h_data_all[tf0] = list(tf_data)
        self.hdr_df = pd.DataFrame(h_data_all)
        # print(f'{len(all_hks)} UNUSED HEADERS: {all_hks}')

    def run(self):
        print("Loading Data...")
        self.load_data()
        print("Loading Header...")
        self.load_header()
        self.signal.finished.emit()
