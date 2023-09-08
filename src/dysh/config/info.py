from getpass import getuser
from pathlib import Path
import os, socket
import platform
import glob

class SystemInfo:

    def __init__(self):
        self.username = self.get_username()
        self.hostname = self.get_hostname()
        self.domain = self.get_domain()
        self.observatory = self.get_observatory()

    def get_username(self):
        return getuser()
    
    def get_hostname(self):
        return socket.gethostname()

    def get_domain(self):
        return socket.gethostbyname_ex(self.hostname)[0]
    
    def get_observatory(self):
        if self.domain.endswith("gb.nrao.edu"):
            observatory = "GBO"
        else:
            observatory = "Unknown"
        return observatory
    
    def get_os(self):
        return platform.uname()[0]

class Paths:
    def __init__(self):
        self.dysh_base = Path(__file__).resolve().parent.parent.parent
        self.data_dirs = []

    def add_data_dir(self, dir):
        if dir not in self.data_dirs:
            self.data_dirs.append(dir)
        else:
            print("Path already in data_dirs")

    def load_gbo(self):
        self.add_data_dir("/home/gbtdata")

