import os
os.environ['QT_API'] = 'pyqt5'
import sip
sip.setapi("QString", 2)
sip.setapi("QVariant", 2)

from PyQt5.QtCore import *
from qtconsole.rich_ipython_widget import RichIPythonWidget
from qtconsole.inprocess import QtInProcessKernelManager
from IPython.lib import guisupport

class QIPythonConsoleWidget(RichIPythonWidget):
    """ Convenience class for a live IPython console widget. We can replace the standard banner using the customBanner argument"""
    def __init__(self, *args, **kwargs):
        super(QIPythonConsoleWidget, self).__init__(*args, **kwargs)

        self.makeBanner()
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel()
        self.kernel_manager.kernel.gui = 'qt4'
        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()
           
        self.exit_requested.connect(self.stop)

    def stop(self):
        self.kernel_client.stop_channels()
        self.kernel_manager.shutdown_kernel()
        guisupport.get_app_qt4().exit() 

    def makeBanner(self):
        self.banner = "Welcome to the Dysh console!\n"

    def pushVariables(self,variableDict):
        """ Given a dictionary containing name / value pairs, push those variables to the IPython console widget """
        self.kernel_manager.kernel.shell.push(variableDict)
    
    def clearTerminal(self):
        """ Clears the terminal """
        self._control.clear()    
    
    def printText(self,text):
        """ Prints some plain text to the console """
        self._append_plain_text(text)        
    
    def executeCommand(self, command):
        """ Execute a command in the frame of the console widget """
        self._execute(command, False)