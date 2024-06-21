# PACKAGE IMPORTS
import sys
import pyqtgraph as pg
from IPython.display import display
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from pyqtgraph import PlotWidget
from pyqtgraph.jupyter import GraphicsLayoutWidget
from qt_material import apply_stylesheet
from screeninfo import get_monitors

# DYSH IMPORTS
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.plot.renderer import Renderer

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

        nomask = self._spectrum.mask
        print("MASK: ", nomask)
        self.plot(sa[~nomask], sf[~nomask], color="white")
        self.plot(sa[nomask], sf[nomask], color="red")

    def get_kwargs(self, **kwargs):
        self.lw =  self._plot_kwargs['linewidth']
        self.xunit = self._plot_kwargs["xaxis_unit"]
        self.yunit = self._plot_kwargs["yaxis_unit"]

class iSingleSpectrum:
    """ Spectrum Plot """
    def __init__(self,spectrum,**kwargs):
        #self.app = pg.mkQApp()
        self.win = GraphicsLayoutWidget(css_width="1000px", css_height="600px")
        self._spectrum = spectrum
        self.add_data()
        self.add_roi()

    def config(self):
        # [TODO] connect spec_num to the hline value
        spec_num = 0
        self.setLabels(left='Intensity', bottom='Frequency')
        self.setTitle(f"Spectrum {spec_num}")

    def update_data(self, spectrum):
        self.clear()
        self._spectrum = spectrum
        self.add_data()

    def add_roi(self):
        "Add a ROI"
        self.roi_rect = pg.RectROI(pos=(1.4E9,-1), size=(5E6,2), pen=pg.mkPen('r', width=2), movable=True, resizable=True)
        self.myplot.addItem(self.roi_rect)

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

class SpectrumROI:

    def __init__(self, roi: pg.RectROI, n: int, xlim: tuple, ylim: tuple):
        self.roi = roi
        self.xlim = xlim
        self.ylim = ylim
        self.roi_pos, self.roi_size = self.lim2pos(xlim, ylim)
        self.make_roi()
        self.make_panel()

    def lim2pos(self, xlim: tuple, ylim: tuple):
        """xy limits to corner position"""
        roi_pos = (xlim[0], ylim[0])
        roi_size = (xlim[1] - xlim[0], ylim[1] - ylim[0])
        return roi_pos, roi_size

    def pos2lim(self, roi_pos: tuple, roi_size: tuple):
        """corner position to xy limits"""
        xlim = (roi_pos[0], roi_pos[0] + roi_size[0])
        ylim = (roi_pos[1], roi_pos[1] + roi_size[1])
        return xlim, ylim

    def make_roi(self):
        # self.roi = pg.RectROI(pos=self.roi_pos, size=self.roi_size, pen=pg.mkPen('r', width=2), movable=True, resizable=True)
        self.roi.setPos(self.roi_pos)
        self.roi.setSize(self.roi_size)
        self.roi.setPen(pg.mkPen('r', width=2))
        self.roi.sigRegionChangeFinished.connect(self.get_ROI)

    def make_panel(self):
        self.roi_xmin_label = QLabel("X Minimum: ")
        self.roi_xmax_label = QLabel("X Maximum: ")
        self.roi_ymin_label = QLabel("Y Minimum: ")
        self.roi_ymax_label = QLabel("Y Maximum: ")

        self.roi_xmin_edit = QLineEdit()
        self.roi_xmax_edit = QLineEdit()
        self.roi_ymin_edit = QLineEdit()
        self.roi_ymax_edit = QLineEdit()

        self.clear_txt_button = QPushButton("Clear Selection")
        self.clear_txt_button.clicked.connect(self.get_ROI)

        self.move_roi_button = QPushButton("Update Selection")
        self.move_roi_button.clicked.connect(self.set_ROI)

        self.select_button = QPushButton("Save Selection")
        self.select_button.clicked.connect(self.get_ROI)

        self.panel_widget = QWidget()
        self.panel_layout = QGridLayout()
        self.panel_widget.setLayout(self.panel_layout)

        self.panel_layout.addWidget(self.roi_xmin_label, 0, 0, 1, 1)
        self.panel_layout.addWidget(self.roi_xmax_label, 1, 0, 1, 1)
        self.panel_layout.addWidget(self.roi_ymin_label, 2, 0, 1, 1)
        self.panel_layout.addWidget(self.roi_ymax_label, 3, 0, 1, 1)

        self.panel_layout.addWidget(self.roi_xmin_edit, 0, 1, 1, 1)
        self.panel_layout.addWidget(self.roi_xmax_edit, 1, 1, 1, 1)
        self.panel_layout.addWidget(self.roi_ymin_edit, 2, 1, 1, 1)
        self.panel_layout.addWidget(self.roi_ymax_edit, 3, 1, 1, 1)

        self.panel_layout.addWidget(self.clear_txt_button, 4, 0, 1, 1)
        self.panel_layout.addWidget(self.move_roi_button, 4, 1, 1, 1)

        self.panel_layout.addWidget(self.select_button, 5, 0, 1, 2)

    def get_ROI(self):
        """ Gets the four corners of the ROI """
        # self.roi_xmin = self.roi.pos()[0]
        # self.roi_xmax = self.roi.pos()[0] + self.roi.size()[0]
        # self.roi_ymin = self.roi.pos()[1]
        # self.roi_ymax = self.roi.pos()[1] + self.roi.size()[1]
        print("Getting ROI")
        self.xlim, self.ylim = self.pos2lim(self.roi.pos(), self.roi.size())
        self.roi_xmin_edit.setText(str(self.xlim[0]))
        self.roi_xmax_edit.setText(str(self.xlim[1]))
        self.roi_ymin_edit.setText(str(self.ylim[0]))
        self.roi_ymax_edit.setText(str(self.ylim[1]))

    def set_ROI(self):
        """ Sets the ROI position based on the text boxes """
        print("Setting ROI")
        self.roi_xmin = float(self.roi_xmin_edit.text())
        self.roi_ymin = float(self.roi_ymin_edit.text())        
        self.roi_xmax = float(self.roi_xmax_edit.text())
        self.roi_ymax = float(self.roi_ymax_edit.text())

        self.roi_pos, self.roi_size = self.lim2pos((self.roi_xmin, self.roi_xmax), (self.roi_ymin, self.roi_ymax))

        self.roi.setPos(self.roi_pos)
        self.roi.setSize(self.roi_size)

class DyshMainWindow(QMainWindow):
    """The main window of the GUI"""

    def __init__(self, spectrum, fpath=None):
        """Initializes the main window"""
        super(DyshMainWindow, self).__init__()

        self.spectrum = spectrum

        self.setWindowTitle("Dysh iPlotter")
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
        """ Create the main panel """
        self._init_UI()

    def _init_UI(self):
        """Creates the skeleton structure of the GUI"""
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)
        self._init_plots()
        self._init_roi_tabs()
        self.main_layout.addWidget(self.spec_plot, 0, 0, 1, 2)

    def _init_roi_tabs(self):
        """ Creates the tabs with the ROI(s) """

        self.roi_tabs = QTabWidget()
        self.add_ROI_tab(1, (1.4E9, 1.405E9), (-1,1))
        self.add_ROI_tab(2, (1.41E9, 1.415E9), (-1,1))

        self.main_layout.addWidget(self.roi_tabs, 0, 2, 1, 1)
        

    def add_ROI_tab(self, n: int, xlim: tuple, ylim: tuple):
        roi_rect_n = pg.RectROI(pos=(0,0), size=(1,1), movable=True, resizable=True)
        self.add_ROI(roi_rect_n, n, xlim, ylim)

    def add_ROI(self, roi: pg.RectROI, n: int, xlim: tuple, ylim: tuple):
        """Adds the ROI(s) to the window"""
        roi_n = SpectrumROI(roi, n, xlim, ylim)
        roi_tab_n = roi_n.panel_widget
        self.roi_tabs.addTab(roi_tab_n, f"ROI {n}")
        self.spec_plot.addItem(roi)

    def add_ROI_nb(self, roi: pg.RectROI, n: int, xlim: tuple, ylim: tuple):
        """Adds the ROI(s) to the window"""
        # roi_n = SpectrumROI(roi, n, xlim, ylim)
        # roi_tab_n = roi_n.panel_widget
        # self.roi_tabs.addTab(roi_tab_n, f"ROI {n}")
        self.spec_plot.addItem(roi)

    def _init_plots(self):
        """Creates the plot canvases"""
        self.spec_plot = SingleSpectrum(self.spectrum)

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

class SpectrumSelectApp(QApplication):
    def __init__(self, spectrum, *args):
        QApplication.__init__(self, *args)
        self.main = DyshMainWindow(spectrum)
        self.main.show()

class SpectrumSelect:
    def __init__(self, spectrum):
        self.spectrum = spectrum
        self._init_outputs()
        self.makePlot()
    
    def _init_outputs(self):
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None

    def makePlot(self):
        rend = Renderer()
        print(rend.info())
        if rend.render_type == "notebook":
            plot = iSingleSpectrum(self.spectrum)
            display(plot.win)
        else:
            app = SpectrumSelectApp(self.spectrum, sys.argv)
            apply_stylesheet(app, theme="dark_purple.xml")
            app.exec_()

if __name__ == "__main__":
    from dysh.fits.gbtfitsload import GBTFITSLoad

    filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
    # filename = "E:/Code/GitHub/Work/forks/dysh/test_file.fits"
    sdfits = GBTFITSLoad(filename)
    psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
    ta = psscan.timeaverage(weights="tsys")
    ta.mask[0:3000] = True
    print("LEN: ", len(ta.mask))

    my_selection = SpectrumSelect(ta)
    print(my_selection.xmin)
