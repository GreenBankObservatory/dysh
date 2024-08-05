import sys
import threading

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig, self.ax = plt.subplots()
        super().__init__(self.fig)
        self.setParent(parent)
        self.plot()

    def plot(self):
        self.ax.clear()
        self.ax.plot([0, 1, 2, 3], [10, 1, 20, 3], "r-")
        self.ax.set_title("Interactive Plot")
        self.draw()

    def update_plot(self, x, y):
        self.ax.clear()
        self.ax.plot(x, y, "r-")
        self.ax.set_title("Updated Plot")
        self.draw()


class Communicator(QObject):
    update_signal = pyqtSignal(list, list)


class MainWindow(QMainWindow):
    def __init__(self, communicator):
        super().__init__()
        self.setWindowTitle("PyQt5 Matplotlib IPython Example")
        self.setGeometry(100, 100, 800, 600)

        self.plot_canvas = PlotCanvas(self)
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
