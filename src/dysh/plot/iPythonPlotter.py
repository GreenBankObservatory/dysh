import sys
import threading

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget

from dysh.plot.canvas import FigureCanvas1D


class Communicator(QObject):
    update_signal = pyqtSignal(list, list)
    signal_update_kwargs = pyqtSignal(dict)


class MainWindow(QMainWindow):
    def __init__(self, communicator):
        super().__init__()
        self.setWindowTitle("PyQt5 Matplotlib IPython Example")
        self.setGeometry(100, 100, 800, 600)

        self.plot_canvas = FigureCanvas1D(self)
        self.communicator = communicator
        self.communicator.update_signal.connect(self.plot_canvas.update_plot)

        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.addWidget(self.plot_canvas)
        self.setCentralWidget(widget)


def start_qt_app(communicator):
    app = QApplication(sys.argv)
    main_window = MainWindow(communicator)
    main_window.show()
    app.exec_()
