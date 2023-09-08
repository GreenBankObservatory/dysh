# PACKAGE IMPORTS
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import sys, os, psutil, getpass, socket, json
import numpy as np
import pyqtgraph as pg
from astropy.io import fits
from time import time
import pandas as pd
import argparse
from screeninfo import get_monitors
from qt_material import apply_stylesheet

class PanelGeneralInfo(QWidget):
    def __init__(self):
        super().__init__()
        self.make_main_panel()

    def make_main_panel(self):
        self.panel_layout = QGridLayout()
        self.setLayout(self.panel_layout)

        self.panel_title = QLabel("General Information")
        self.panel_layout.addWidget(self.panel_title, 0, 0, 1, 2)

        self.make_radio_derive()
        self.panel_layout.addWidget(QLabel("Derive"), 1, 0, 1, 1)
        self.panel_layout.addWidget(self.radio_set_derive, 1, 1, 1, 1)

        self.make_radio_sensitivity_units()
        self.panel_layout.addWidget(QLabel("Sensitivity Units"), 2, 0, 1, 1)
        self.panel_layout.addWidget(self.radio_set_sensitivity_units, 2, 1, 1, 1)

        self.make_text_edit_sensitivity()
        self.panel_layout.addWidget(QLabel("Desired Sensitivity (1-sigma)"), 3, 0, 1, 1)
        self.panel_layout.addWidget(self.text_edit_sensitivity, 3, 1, 1, 1)

    def make_radio_derive(self):
        self.radio_set_derive = QWidget()
        self.radio_set_derive_layout = QVBoxLayout()
        self.radio_set_derive.setLayout(self.radio_set_derive_layout)

        self.rs_derive_0 = QRadioButton("Observing Time from Desired Sensitivity")
        self.rs_derive_1 = QRadioButton("Sensitivity from Observing Time")
        
        self.rs_derive_0.setChecked(True)
        self.rs_derive_1.setChecked(False)
        
        self.rs_derive_0.toggled.connect(lambda:self.rs_derive_btnstate(self.rs_derive_0))
        self.rs_derive_1.toggled.connect(lambda:self.rs_derive_btnstate(self.rs_derive_1))
        
        self.radio_set_derive_layout.addWidget(self.rs_derive_0)
        self.radio_set_derive_layout.addWidget(self.rs_derive_1)

    def rs_derive_btnstate(self,b):
	
        if b.text() == "Observing Time from Desired Sensitivity":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass
                
        if b.text() == "Sensitivity from Observing Time":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

    def make_radio_sensitivity_units(self):
        self.radio_set_sensitivity_units = QWidget()
        self.radio_set_sensitivity_units_layout = QVBoxLayout()
        self.radio_set_sensitivity_units.setLayout(self.radio_set_sensitivity_units_layout)

        self.rs_sensitivity_units_0 = QRadioButton("Flux Density (mJy)")
        self.rs_sensitivity_units_1 = QRadioButton("Antenna Temp., Ta (mK)")
        self.rs_sensitivity_units_2 = QRadioButton("Main Beam Temp., Tmb (mK)")
        self.rs_sensitivity_units_3 = QRadioButton("Forward Antenna Temperature, Ta* (mK)")
        self.rs_sensitivity_units_4 = QRadioButton("Radiation Temp., Tr (mK)")
        
        self.rs_sensitivity_units_0.setChecked(True)
        self.rs_sensitivity_units_1.setChecked(False)
        self.rs_sensitivity_units_2.setChecked(False)
        self.rs_sensitivity_units_3.setChecked(False)
        self.rs_sensitivity_units_4.setChecked(False)
        
        self.rs_sensitivity_units_0.toggled.connect(lambda:self.rs_sensitivity_units_btnstate(self.rs_sensitivity_units_0))
        self.rs_sensitivity_units_1.toggled.connect(lambda:self.rs_sensitivity_units_btnstate(self.rs_sensitivity_units_1))
        self.rs_sensitivity_units_2.toggled.connect(lambda:self.rs_sensitivity_units_btnstate(self.rs_sensitivity_units_2))
        self.rs_sensitivity_units_3.toggled.connect(lambda:self.rs_sensitivity_units_btnstate(self.rs_sensitivity_units_3))
        self.rs_sensitivity_units_4.toggled.connect(lambda:self.rs_sensitivity_units_btnstate(self.rs_sensitivity_units_4))
        
        self.radio_set_sensitivity_units_layout.addWidget(self.rs_sensitivity_units_0)
        self.radio_set_sensitivity_units_layout.addWidget(self.rs_sensitivity_units_1)
        self.radio_set_sensitivity_units_layout.addWidget(self.rs_sensitivity_units_2)
        self.radio_set_sensitivity_units_layout.addWidget(self.rs_sensitivity_units_3)
        self.radio_set_sensitivity_units_layout.addWidget(self.rs_sensitivity_units_4)

    def rs_sensitivity_units_btnstate(self,b):
	
        if b.text() == "Flux Density (mJy)":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass
				
        if b.text() == "Antenna Temp., Ta (mK)":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

        if b.text() == "Main Beam Temp., Tmb (mK)":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

        if b.text() == "Forward Antenna Temperature, Ta* (mK)":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

        if b.text() == "Radiation Temp., Tr (mK)":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

    def make_text_edit_sensitivity(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_sensitivity = QLineEdit()
        self.text_edit_sensitivity.setValidator(QIntValidator())
        self.text_edit_sensitivity.textChanged.connect(self.validate_text_edit_sensitivity)

    def validate_text_edit_sensitivity(self):
        """ Validates the input in self.text_edit_sensitivity """
        pass

class PanelHardwareInfo(QWidget):
    def __init__(self):
        super().__init__()
        self.make_main_panel()

    def make_main_panel(self):
        self.panel_hi_layout = QGridLayout()
        self.setLayout(self.panel_hi_layout)

        self.panel_hi_title = QLabel("Hardware Information")
        self.panel_hi_layout.addWidget(self.panel_hi_title)
        
        self.panel_hi_info = QLabel(
            """
            Answer questions from top to bottom. 
            If you change a question that was answered previously, 
            check all answers that follow. 
            Some answers will dictate the answer for other questions.
            """
            )
        
        self.panel_hi_layout.addWidget(self.panel_hi_info, 1, 0, 1, 1)

        self.make_combo_backend()
        self.panel_hi_layout.addWidget(QLabel("Backend"), 2, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_backend, 2, 1, 1, 1)

        self.make_combo_mode()
        self.panel_hi_layout.addWidget(QLabel("Mode"), 3, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_mode, 3, 1, 1, 1)

        self.make_combo_rx()
        self.panel_hi_layout.addWidget(QLabel("Receiver"), 4, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_rx, 4, 1, 1, 1)
        self.panel_hi_layout.addWidget(QLabel("Beams"), 5, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_beams, 5, 1, 1, 1)
        self.panel_hi_layout.addWidget(QLabel("Polarization"), 6, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_pol, 6, 1, 1, 1)
        self.panel_hi_layout.addWidget(QLabel("Bandwidth per Window (MHz)"), 7, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_bpw, 7, 1, 1, 1)
        self.panel_hi_layout.addWidget(QLabel("Number of Spectral Windows"), 8, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_n_spec, 8, 1, 1, 1)
        self.panel_hi_layout.addWidget(QLabel("Switching Mode"), 9, 0, 1, 1)
        self.panel_hi_layout.addWidget(self.combo_sw_mode, 9, 1, 1, 1)

    def make_combo_backend(self):
        self.combo_backend = QComboBox()
        self.combo_backend.addItem("VErsatile GB Astronomical Spectrometer")
        self.combo_backend.setEnabled(False)

    def make_combo_mode(self):
        self.combo_mode = QComboBox()
        self.combo_mode.addItem("Spectral Line")
        self.combo_mode.setEnabled(False)

    def make_combo_rx(self):
        f = open("/Users/victoriacatlett/Documents/Code/repos/github/dysh/gui/static/gbt/gbt-rx.json")
        self.rx_dict = json.load(f)
        
        self.combo_rx = QComboBox()
        for k in self.rx_dict:
            self.combo_rx.addItem(k)

        self.make_combo_beams()
        self.make_combo_pol()
        self.make_combo_bpw()
        self.make_combo_n_spec()
        self.make_combo_sw_mode()

        self.combo_rx.currentIndexChanged.connect(self.update_combo_rx)

        self.update_combo_rx()

    def make_combo_beams(self):
        self.combo_beams = QComboBox()

    def make_combo_pol(self):
        self.combo_pol = QComboBox()
    
    def make_combo_bpw(self):
        self.combo_bpw = QComboBox()

    def make_combo_n_spec(self):
        self.combo_n_spec = QComboBox()

    def make_combo_sw_mode(self):
        self.combo_sw_mode = QComboBox()

    def update_combo_rx(self):
        # [TODO] Load the RX info from a JSON file
        self.update_combo_beams()
        self.update_combo_pol()
        self.update_combo_bpw()
        self.update_combo_n_spec()
        self.update_combo_sw_mode()
    
    def update_combo_beams(self):
        self.combo_beams.clear()
        self.combo_beams.setEnabled(True)
        ct = self.combo_rx.currentText()
        for bi in self.rx_dict[ct]["beams"]:
            self.combo_beams.addItem(str(bi))
        if len(self.rx_dict[ct]["beams"]) < 2:
            self.combo_beams.setEnabled(False)

    def update_combo_pol(self):
        self.combo_pol.clear()
        self.combo_pol.setEnabled(True)
        ct = self.combo_rx.currentText()
        for pi in self.rx_dict[ct]["pol"]:
            self.combo_pol.addItem(str(pi))
        if len(self.rx_dict[ct]["pol"]) < 2:
            self.combo_pol.setEnabled(False)
    
    def update_combo_bpw(self):
        self.combo_bpw.clear()
        self.combo_bpw.setEnabled(True)
        ct = self.combo_rx.currentText()
        for bi in self.rx_dict[ct]["BpW"]:
            self.combo_bpw.addItem(str(bi))
        if len(self.rx_dict[ct]["BpW"]) < 2:
            self.combo_bpw.setEnabled(False)
    
    def update_combo_n_spec(self):
        self.combo_n_spec.clear()
        self.combo_n_spec.setEnabled(True)
        ct = self.combo_rx.currentText()
        for ni in self.rx_dict[ct]["n_spec"]:
            self.combo_n_spec.addItem(str(ni))
        if len(self.rx_dict[ct]["n_spec"]) < 2:
            self.combo_n_spec.setEnabled(False)
    
    def update_combo_sw_mode(self):
        # [TODO] Add logic for reducing number for ARGUS with 4, 8, or 16 beams
        self.combo_sw_mode.clear()
        self.combo_sw_mode.setEnabled(True)
        ct = self.combo_rx.currentText()
        for si in self.rx_dict[ct]["sw_mode"]:
            self.combo_sw_mode.addItem(str(si))
        if len(self.rx_dict[ct]["sw_mode"]) < 2:
            self.combo_sw_mode.setEnabled(False)

class PanelSourceInfo(QWidget):
    def __init__(self):
        super().__init__()
        self.make_main_panel()

    def make_main_panel(self):
        self.panel_layout = QGridLayout()
        self.setLayout(self.panel_layout)

        self.panel_title = QLabel("Source Information")
        self.panel_layout.addWidget(self.panel_title, 0, 0, 1, 2)

        self.make_radio_freq_spec()
        self.panel_layout.addWidget(QLabel("Frequency Specified in the:"), 1, 0, 1, 2)
        self.panel_layout.addWidget(self.radio_set_freq_spec, 1, 1, 1, 2)

        self.make_text_edit_freq()
        self.panel_layout.addWidget(QLabel("Rest Frequency (MHz):"), 2, 0, 1, 2)
        self.panel_layout.addWidget(self.text_edit_freq, 2, 1, 1, 2)

        self.make_combo_doppler_correction()
        self.panel_layout.addWidget(QLabel("Doppler Correction:"), 3, 0, 1, 2)
        self.panel_layout.addWidget(self.combo_doppler_correction, 3, 1, 1, 2)

        self.make_text_edit_souvel()
        self.panel_layout.addWidget(QLabel("Source Velocity (km/s):"), 4, 0, 1, 2)
        self.panel_layout.addWidget(self.text_edit_souvel, 4, 1, 1, 2)

        self.make_text_edit_soudiam()
        self.panel_layout.addWidget(QLabel("Source Diameter (arc minutes):"), 5, 0, 1, 2)
        self.panel_layout.addWidget(self.text_edit_soudiam, 5, 1, 1, 2)

        self.panel_layout.addWidget(QLabel("Source Contribution Corrections"), 6, 0, 1, 1)

        self.make_radio_contrib()
        self.panel_layout.addWidget(QLabel("Source Contribution to System Temperature:"), 7, 0, 1, 1)
        self.panel_layout.addWidget(self.radio_set_contrib, 7, 1, 1, 1)

        self.make_text_edit_contrib_k()
        self.panel_layout.addWidget(QLabel("Contribution (K):"), 8, 0, 1, 1)
        self.panel_layout.addWidget(self.text_edit_contrib_k, 8, 1, 1, 1)

        self.make_text_edit_soudec()
        self.panel_layout.addWidget(QLabel("Source Declination (Deg):"), 9, 0, 1, 2)
        self.panel_layout.addWidget(self.text_edit_soudec, 9, 1, 1, 2)

        self.make_text_edit_min_el()
        self.panel_layout.addWidget(QLabel("Minimum Elevation (Deg):"), 10, 0, 1, 2)
        self.panel_layout.addWidget(self.text_edit_min_el, 10, 1, 1, 2)

    def make_radio_freq_spec(self):
        self.radio_set_freq_spec = QWidget()
        self.radio_set_freq_spec_layout = QVBoxLayout()
        self.radio_set_freq_spec.setLayout(self.radio_set_freq_spec_layout)

        self.rs_freq_spec_0 = QRadioButton("Topocentric Frame")
        self.rs_freq_spec_1 = QRadioButton("Rest Frame")
        
        self.rs_freq_spec_0.setChecked(True)
        self.rs_freq_spec_1.setChecked(False)
        
        self.rs_freq_spec_0.toggled.connect(lambda:self.rs_freq_spec_btnstate(self.rs_freq_spec_0))
        self.rs_freq_spec_1.toggled.connect(lambda:self.rs_freq_spec_btnstate(self.rs_freq_spec_1))
        
        self.radio_set_freq_spec_layout.addWidget(self.rs_freq_spec_0)
        self.radio_set_freq_spec_layout.addWidget(self.rs_freq_spec_1)

    def rs_freq_spec_btnstate(self,b):
	
        if b.text() == "Topocentric Frame":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass
                
        if b.text() == "Rest Frame":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

    def make_text_edit_freq(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_freq = QLineEdit()
        self.text_edit_freq.setValidator(QDoubleValidator())
        self.text_edit_freq.textChanged.connect(self.validate_text_edit_freq)

    def validate_text_edit_freq(self):
        """ Validates the input in self.text_edit_freq """
        pass

    def make_combo_doppler_correction(self):
        self.combo_doppler_correction = QComboBox()
        self.combo_doppler_correction.currentIndexChanged.connect(self.update_combo_doppler_correction)

        self.combo_doppler_correction.addItem("Radio")
        self.combo_doppler_correction.addItem("Optical")
        self.combo_doppler_correction.addItem("Redshift")

    def update_combo_doppler_correction(self):
        ct = self.combo_doppler_correction.currentText()
        print(ct)

    def make_text_edit_souvel(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_souvel = QLineEdit()
        self.text_edit_souvel.setValidator(QDoubleValidator())
        self.text_edit_souvel.textChanged.connect(self.validate_text_edit_souvel)

    def validate_text_edit_souvel(self):
        """ Validates the input in self.text_edit_souvel """
        pass

    def make_text_edit_soudiam(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_soudiam = QLineEdit()
        self.text_edit_soudiam.setValidator(QDoubleValidator())
        self.text_edit_soudiam.textChanged.connect(self.validate_text_edit_soudiam)

    def validate_text_edit_soudiam(self):
        """ Validates the input in self.text_edit_soudiam """
        pass

    def make_radio_contrib(self):
        self.radio_set_contrib = QWidget()
        self.radio_set_contrib_layout = QVBoxLayout()
        self.radio_set_contrib.setLayout(self.radio_set_contrib_layout)

        self.rs_contrib_0 = QRadioButton("User Estimated Correction")
        self.rs_contrib_1 = QRadioButton("Internal Galactic Model")
        
        self.rs_contrib_0.setChecked(True)
        self.rs_contrib_1.setChecked(False)
        
        self.rs_contrib_0.toggled.connect(lambda:self.rs_contrib_btnstate(self.rs_contrib_0))
        self.rs_contrib_1.toggled.connect(lambda:self.rs_contrib_btnstate(self.rs_contrib_1))
        
        self.radio_set_contrib_layout.addWidget(self.rs_contrib_0)
        self.radio_set_contrib_layout.addWidget(self.rs_contrib_1)

    def rs_contrib_btnstate(self,b):
	
        if b.text() == "User Estimated Correction":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass
                
        if b.text() == "Internal Galactic Model":
            if b.isChecked() == True:
                print(b.text()+" is selected")
            else:
                pass

    def make_text_edit_contrib_k(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_contrib_k = QLineEdit()
        self.text_edit_contrib_k.setValidator(QDoubleValidator())
        self.text_edit_contrib_k.textChanged.connect(self.validate_text_edit_contrib_k)

    def validate_text_edit_contrib_k(self):
        """ Validates the input in self.text_edit_contrib_k """
        pass

    def make_text_edit_soudec(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_soudec = QLineEdit()
        self.text_edit_soudec.setValidator(QDoubleValidator())
        self.text_edit_soudec.textChanged.connect(self.validate_text_edit_soudec)

    def validate_text_edit_soudec(self):
        """ Validates the input in self.text_edit_soudec """
        pass

    def make_text_edit_min_el(self):
        """ Text edit to type the desired sensitivity """
        self.text_edit_min_el = QLineEdit()
        self.text_edit_min_el.setValidator(QDoubleValidator())
        self.text_edit_min_el.textChanged.connect(self.validate_text_edit_min_el)

    def validate_text_edit_min_el(self):
        """ Validates the input in self.text_edit_min_el """
        pass

class PanelDataReduction(QWidget):
    def __init__(self):
        super().__init__()
        self.make_main_panel()

    def make_main_panel(self):
        self.panel_layout = QGridLayout()
        self.setLayout(self.panel_layout)

        self.panel_title = QLabel("Data Reduction")
        self.panel_layout.addWidget(self.panel_title, 0, 0, 1, 2)

class SensitivityCalcForm(QMainWindow):
    """ The main window of the GUI """
    def __init__(self):
        """ Initializes the main window """
        super(SensitivityCalcForm, self).__init__()
        self.setWindowTitle("Sensitivity Calculator GUI")

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

        self.make_panel_general_info()
        self.make_panel_hardware_info()
        self.make_panel_source_info()
        self.make_panel_data_reduction()

        self.main_layout.addWidget(self.panel_general_info, 0, 0, 1, 1)
        self.main_layout.addWidget(self.panel_hardware_info, 1, 0, 1, 1)
        self.main_layout.addWidget(self.panel_source_info, 0, 1, 1, 1)
        self.main_layout.addWidget(self.panel_data_reduction, 1, 1, 1, 1)

    def make_panel_general_info(self):
        self.panel_general_info = PanelGeneralInfo()

    def make_panel_hardware_info(self):
        self.panel_hardware_info = PanelHardwareInfo()

    def make_panel_source_info(self):
        self.panel_source_info = PanelSourceInfo()

    def make_panel_data_reduction(self):
        self.panel_data_reduction = PanelDataReduction()

def start():
    app = QApplication(sys.argv)
    apply_stylesheet(app, theme='dark_purple.xml')
    win = SensitivityCalcForm()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()