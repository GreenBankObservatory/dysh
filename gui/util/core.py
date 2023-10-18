from PyQt5.QtCore import QObject, QThread, pyqtSignal
from threading import Thread
import sys, os

class ThreadCallbacks:
    def progress(future):
        pass
        # print('.', end='', flush=True)

class DyshWorker(QObject, Thread):
    """Thread to run a function that returns a value"""

    def __init__(self, group=None, target=None, name=None, args=(), kwargs={}):
        self.finished = pyqtSignal()
        self.progress = pyqtSignal(int)
        QObject.__init__(self)
        Thread.__init__(self, group=group, target=target, name=name, args=args, kwargs=kwargs)
        self._return = None
        # print(self.__dir__())

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self):
        Thread.join(self)
        return self._return