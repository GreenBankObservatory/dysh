import sys

import pyqtgraph as pg
from IPython.display import display
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from pyqtgraph import PlotWidget
from qt_material import apply_stylesheet
from screeninfo import get_monitors

# DYSH IMPORTS
from dysh.plot.renderer import Renderer

# import threading
# import time


class SingleSpectrum(PlotWidget):
    """Spectrum Plot"""

    def __init__(self, spectrum):
        """Init the spectrum plot"""
        super().__init__()
        self._spectrum = spectrum
        self.reset_kwargs()
        self.config()
        self.add_data()

    def config(self):
        """Get plot configuration"""
        self.setLabels(
            left=f'{self._plot_kwargs["ylabel"]} ({self._plot_kwargs["yaxis_unit"]})',
            bottom=f'{self._plot_kwargs["xlabel"]} ({self._plot_kwargs["xaxis_unit"]})',
        )
        self.setTitle(self._plot_kwargs["title"])

    def update_data(self, spectrum):
        """Update data if spectrum changes"""
        self.clear()
        self._spectrum = spectrum
        self.add_data()

    def add_data(self):
        """Add the data to the plot"""

        s = self._spectrum
        sa = s.spectral_axis
        sf = s.flux

        nomask = self._spectrum.mask == False
        self.plot(sa[~nomask], sf[~nomask], color=self._plot_kwargs["line_color"])
        self.plot(sa[nomask], sf[nomask], color=self._plot_kwargs["line_color_masked"])

    def set_kwargs(self, key, val):
        """Get plot kwargs"""
        self._plot_kwargs[key] = val

    def reset_kwargs(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = {
            "xmin": None,
            "xmax": None,
            "ymin": None,
            "ymax": None,
            "xlabel": "Frequency",
            "ylabel": "Flux",
            "xaxis_unit": str(self._spectrum.spectral_axis.unit),
            "yaxis_unit": str(self._spectrum.unit),
            "grid": False,
            "figsize": None,
            "xlabel_color": "white",
            "ylabel_color": "white",
            "tick_color": "white",
            "spine_color": "white",
            "linewidth": 2.0,
            "linestyle": "steps-mid",
            "markersize": 8,
            "line_color": "white",
            "line_color_masked": "red",
            "title": "Dysh plot",
            "title_color": "white",
            "aspect": "auto",
            "bbox_to_anchor": None,
            "loc": "best",
            "legend": None,
            "show_baseline": True,
            "test": False,
            "titlecolor": "white",
            "facecolor": "black",
            "edgecolor": "none",
        }


class DyshMainWindow(QMainWindow):
    """The main window of the GUI"""

    def __init__(self, spectrum, add_roi=True):
        """Initializes the main window"""
        super().__init__()

        self.spectrum = spectrum
        self.bool_roi = add_roi

        self.setWindowTitle("Dysh iPlotter")
        self._init_geometry(0.8)

        self.info_threads()
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
                self.xpos = int(m.x + (m.width * (1 - mult)) / 2)
                self.ypos = int(m.y + (m.height * (1 - mult)) / 2)
        self.setGeometry(self.xpos, self.ypos, self.width, self.height)

    def info_threads(self):
        """Updates info on available threads"""
        self.threadCountActive = QThreadPool.globalInstance().activeThreadCount()
        self.threadCountTotal = QThreadPool.globalInstance().maxThreadCount()
        # print(f"You are using {self.threadCountActive} of the {self.threadCountTotal} available QThreads")

    def _init_UI(self):
        """Creates the skeleton structure of the GUI"""
        self.main_widget = QWidget()
        self.main_layout = QGridLayout()
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)
        self._init_plots()
        self.main_layout.addWidget(self.spec_plot, 0, 0, 1, 1)

        if self.bool_roi:
            self._init_roi_panel()
            self.main_layout.addWidget(self.roi_panel_widget, 1, 0, 1, 1)

    def _init_plots(self):
        """Creates the plot canvases"""
        self.spec_plot = SingleSpectrum(self.spectrum)
        if self.bool_roi:
            self._init_roi()

    def _init_roi_panel(self):
        self.roi_panel_widget = QWidget()
        self.roi_panel_layout = QGridLayout()
        self.roi_panel_widget.setLayout(self.roi_panel_layout)

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

        # self.select_button = QPushButton("Save Selection")
        # self.select_button.clicked.connect(self.get_ROI)

        self.roi_panel_layout.addWidget(self.roi_xmin_label, 0, 0, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_xmax_label, 1, 0, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_ymin_label, 2, 0, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_ymax_label, 3, 0, 1, 1)

        self.roi_panel_layout.addWidget(self.roi_xmin_edit, 0, 1, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_xmax_edit, 1, 1, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_ymin_edit, 2, 1, 1, 1)
        self.roi_panel_layout.addWidget(self.roi_ymax_edit, 3, 1, 1, 1)

        self.roi_panel_layout.addWidget(self.clear_txt_button, 4, 0, 1, 1)
        self.roi_panel_layout.addWidget(self.move_roi_button, 4, 1, 1, 1)

        # self.roi_panel_layout.addWidget(self.select_button, 5, 0, 1, 2)

    def _init_roi(self):
        """Make the ROI"""
        self.roi = pg.RectROI((1.39e9, -1), (1e6, 2), pen=pg.mkPen("r", width=2), movable=True, resizable=True)
        self.spec_plot.addItem(self.roi)
        self.roi.sigRegionChangeFinished.connect(self.get_ROI)
        self.xlim, self.ylim = self.roi_pos2lim(self.roi.pos(), self.roi.size())
        self.roi_xmin = self.xlim[0]
        self.roi_ymin = self.ylim[0]
        self.roi_xmax = self.xlim[1]
        self.roi_ymax = self.ylim[1]

    def get_ROI(self):
        """Gets the four corners of the ROI"""
        self.xlim, self.ylim = self.roi_pos2lim(self.roi.pos(), self.roi.size())
        self.roi_xmin_edit.setText(str(self.xlim[0]))
        self.roi_xmax_edit.setText(str(self.xlim[1]))
        self.roi_ymin_edit.setText(str(self.ylim[0]))
        self.roi_ymax_edit.setText(str(self.ylim[1]))

    def set_ROI(self):
        """Sets the ROI position based on the text boxes"""
        self.roi_xmin = float(self.roi_xmin_edit.text())
        self.roi_ymin = float(self.roi_ymin_edit.text())
        self.roi_xmax = float(self.roi_xmax_edit.text())
        self.roi_ymax = float(self.roi_ymax_edit.text())

        self.roi_pos, self.roi_size = self.roi_lim2pos((self.roi_xmin, self.roi_xmax), (self.roi_ymin, self.roi_ymax))

        self.roi.setPos(self.roi_pos)
        self.roi.setSize(self.roi_size)

    def roi_lim2pos(self, xlim: tuple, ylim: tuple):
        """ROI xy limits to corner position"""
        roi_pos = (xlim[0], ylim[0])
        roi_size = (xlim[1] - xlim[0], ylim[1] - ylim[0])
        return roi_pos, roi_size

    def roi_pos2lim(self, roi_pos: tuple, roi_size: tuple):
        """ROI corner position to xy limits"""
        xlim = (roi_pos[0], roi_pos[0] + roi_size[0])
        ylim = (roi_pos[1], roi_pos[1] + roi_size[1])
        return xlim, ylim

    def _clear_all(self):
        """Clear all widgets in the window"""
        while self.main_layout.count():
            child = self.main_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    # def closeEvent(self, *args):
    #     super(QMainWindow, self).closeEvent(*args)


class SpectrumSelectApp(QApplication):
    def __init__(self, spectrum, *args, **kwargs):
        # [TODO] Handle incoming kwargs
        QApplication.__init__(self, *args)
        self.main = DyshMainWindow(spectrum)


class SpectrumPlot:
    """The spectrum selection"""

    def __init__(self, spectrum, **kwargs):
        """Initialize the spectrum seelctor"""
        # [TODO] Handle incoming kwargs
        self.spectrum = spectrum
        self._init_outputs()

        # Attempt at threading that ALMOST works
        # t = threading.Thread(target=self.makePlot)
        # t.daemon = True
        # t.start()

        # if timeout is not None:
        #     print(f"Timing out after {timeout} second(s)")
        #     time.sleep(timeout)
        #     print("Closing main window")
        #     self.app.main.close()

        # Delete once threading works
        self.makePlot()

    def _init_outputs(self):
        """Initialize the selection object"""
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None

    def makePlot(self):
        """Make the plot based on the renderer"""

        rend = Renderer()

        if rend.render_type == "notebook":
            plot = SingleSpectrum(self.spectrum)
            display(plot.win)
        else:
            print("app started")
            self.app = SpectrumSelectApp(self.spectrum, sys.argv)
            apply_stylesheet(self.app, theme="dark_purple.xml")
            self.app.exec()
            self.xmin = self.app.main.roi_xmin
            self.xmax = self.app.main.roi_xmax
            self.ymin = self.app.main.roi_ymin
            self.ymax = self.app.main.roi_ymax

    def get_selection(self):
        """Return the selection as a dictionary"""
        xunit = self.spectrum.spectral_axis.unit
        yunit = self.spectrum.unit
        self.selection = {
            "x": (self.xmin * xunit, self.xmax * xunit),
            "y": (self.ymin * yunit, self.ymax * yunit),
        }
        return self.selection

    def close(self):
        self.app.quit()
