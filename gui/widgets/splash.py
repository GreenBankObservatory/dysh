import os
import sys
from pathlib import Path

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QMovie
from PyQt5.QtWidgets import QApplication, QMainWindow, QSplashScreen

GUI_BASE_DIR = Path(__file__).resolve().parent.parent
# https://stackoverflow.com/questions/71627508/pyqt-show-animated-gif-while-other-operations-are-running


class Worker(QThread):
    progressChanged = pyqtSignal(int)

    def run(self):
        for count in range(6):
            self.progressChanged.emit(count)
            self.sleep(1)
        self.progressChanged.emit(-1)


class Window(QMainWindow):
    def __init__(self):
        super().__init__()


class SplashScreen(QSplashScreen):
    def __init__(self, flags=0):
        super().__init__(flags=Qt.WindowFlags(flags))
        self.load_gif_dir = os.path.join(GUI_BASE_DIR, "static/img/loading.gif")
        self.movie = QMovie(self.load_gif_dir, parent=self)
        self.movie.frameChanged.connect(self.handleFrameChange)
        self.movie.start()

    def updateProgress(self, count=0):
        if count == 0:
            message = "Starting..."
        elif count > 0:
            message = f"Processing... {count}"
        else:
            message = "Finished!"
        self.showMessage(message, Qt.AlignHCenter | Qt.AlignBottom, Qt.white)

    def handleFrameChange(self):
        pixmap = self.movie.currentPixmap()
        self.setPixmap(pixmap)
        self.setMask(pixmap.mask())


if __name__ == "__main__":
    print("Splash screen!")
    dir_path = os.path.dirname(os.path.realpath(__file__))

    app = QApplication(sys.argv)
    window = Window()

    splash = SplashScreen(Qt.WindowStaysOnTopHint)
    worker = Worker()
    worker.progressChanged.connect(splash.updateProgress)
    worker.finished.connect(lambda: (splash.finish(window), window.show()))
    splash.show()
    worker.start()
    app.exec_()
