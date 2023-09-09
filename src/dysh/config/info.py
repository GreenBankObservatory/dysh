from getpass import getuser
import os, socket
import platform
import glob


class SystemInfo:
    def __init__(self):
        self.USERNAME = self.get_username()
        self.HOSTNAME = self.get_hostname()
        self.DOMAIN = self.get_domain()
        self.OBSERVATORY = self.get_observatory()

    def get_username(self):
        return getuser()

    def get_hostname(self):
        return socket.gethostname()

    def get_domain(self):
        return socket.gethostbyname_ex(self.HOSTNAME)[0]

    def get_observatory(self):
        if self.DOMAIN.endswith("gb.nrao.edu"):
            observatory = "GBO"
        else:
            observatory = "Unknown"
        return observatory

    def get_os(self):
        return platform.uname()[0]
