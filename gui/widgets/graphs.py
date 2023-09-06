from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from pyqtgraph import GraphicsLayoutWidget, PlotWidget, ImageView, InfiniteLine
import sys, os, psutil, getpass, socket

class SpectrumSelectLine(InfiniteLine):
    """ Horizontal Line to select spectrum in waterfall """
    def __init__(self):
        super().__init__(pos=0, angle=0, movable=True)
        self.config()
    
    def config(self):
        self.setPen(color='r', width=2)
        self.setHoverPen(color='r', width=2)

    def connect(self, func):
        self.sigDragged.connect(func)

class WaterfallSpectrum(ImageView):
    """ Waterfall Plot """
    def __init__(self, data):
        super().__init__()
        self.add_data(data)
        self.config()

    def add_data(self, data):
        self.setImage(self.data.T)

    def config(self):
        #self.setLimits(xMin=0, xMax=self.xlim0, yMin=self.ylim0, yMax=0)
        self.setTitle('Waterfall Plot')
        self.setLabels(left='Spectrum', bottom='Frequency')
        self.setColorMap('viridis')

    def add_hline(self, change_func):
        self.hline = SpectrumSelectLine()
        self.hline.connect(change_func)
        self.addItem(self.hline)

class SingleSpectrum(PlotWidget):
    """ Spectrum Plot """
    def __init__(self):
        super().__init__()

    def config(self):
        # [TODO] connect spec_num to the hline value
        spec_num = 0
        self.setLabels(left='Intensity', bottom='Frequency')
        self.setTitle(f"Spectrum {spec_num}")

    def update_data(self, data):
        self.clear()
        self.add_data(data)

    def add_data(self, data):
        self.plot(self.data.x, self.data.y)

class MemoryUsage(PlotWidget):
    """ Memory Usage Plot """
    def __init__(self):
        super().__init__()
        self.username = getpass.getuser()
        self.hostname = socket.gethostname()

    def _init_mem_use_plot(self):
        self.list_mem_time = []
        self.list_pmem_val = []
        self.list_amem_val = []
        self.plotWidget_mem = pg.GraphicsLayoutWidget()
        self.ax_mem = self.plotWidget_mem.addPlot(title=f"RAM Usage ({self.username}@{self.hostname})", row=1, col=0)
        self.ax_mem.setLabels(left='Resident Set Size (MiB)', bottom='Time')
        #self.ax_mem.setLogMode(False, True)

        self.timer_mem = QTimer()
        self.timer_mem.timeout.connect(self.update_mem_use_plot)
        self.timer_mem.start(1000)

    def update_mem_use_plot(self):
        # [TODO] Figure out what to do if this is left open for a long time. The memory lists will get huge. 
        pmem, amem = self.get_memory_usage()
        self.list_pmem_val.append(pmem)
        self.list_amem_val.append(amem)
        self.list_mem_time.append(time())
        self.ax_mem.clear()
        #self.ax_mem.addLegend()
        self.hline_fsize = pg.InfiniteLine(pos=self.fsize, angle=0, movable=False, label="SDFITS File Size", labelOpts={'color':'r'})
        self.hline_fsize.setPen(color="r", width=1)
        self.ax_mem.addItem(self.hline_fsize, name="SDFITS File Size")
        #self.ax_mem.plot(self.list_mem_time, self.list_amem_val, name=f"Available RAM on {self.hostname}")
        self.ax_mem.plot(self.list_mem_time, self.list_pmem_val, name="Current RAM Usage")

    def get_memory_usage(self):
        # [TODO] Verify that this is not underestimating usage bc of caches
        process = psutil.Process(os.getpid())
        pmem = process.memory_info().rss / 1048576 # Bytes to MiB
        amem = psutil.virtual_memory().available / 1048576 # Bytes to MiB
        return pmem, amem