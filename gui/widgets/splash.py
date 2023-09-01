import sys, os
from PyQt5.QtCore import pyqtSignal, Qt, QThread
from PyQt5.QtGui import QMovie
from PyQt5.QtWidgets import QApplication, QSplashScreen, QMainWindow
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
    def __init__(self, filepath, flags=0):
        super().__init__(flags=Qt.WindowFlags(flags))
        self.movie = QMovie(filepath, parent=self)
        self.movie.frameChanged.connect(self.handleFrameChange)
        self.movie.start()

    def updateProgress(self, count=0):
        if count == 0:
            message = 'Starting...'
        elif count > 0:
            message = f'Processing... {count}'
        else:
            message = 'Finished!'
        self.showMessage(
            message, Qt.AlignHCenter | Qt.AlignBottom, Qt.white)
        
    def handleFrameChange(self):
        pixmap = self.movie.currentPixmap()
        self.setPixmap(pixmap)
        self.setMask(pixmap.mask())

if __name__ == "__main__":
    print("Splash screen!")
    dir_path = os.path.dirname(os.path.realpath(__file__))

    app = QApplication(sys.argv)
    window = Window()
    load_path = dir_path + "/img/loading.gif"

    splash = SplashScreen(load_path, Qt.WindowStaysOnTopHint)
    worker = Worker()
    worker.progressChanged.connect(splash.updateProgress)
    worker.finished.connect(
        lambda: (splash.finish(window), window.show()))
    splash.show()
    worker.start()
    app.exec_()