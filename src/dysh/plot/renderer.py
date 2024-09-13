import sys
from pathlib import Path

from IPython import get_ipython

from dysh import __version__ as dysh_version


class EnvironmentInfo:
    def __init__(self):
        self.environment = self.detect_environment()

    def detect_environment(self):
        ip = get_ipython()
        if ip is None:
            return "Python Script"
        else:
            if ip.has_trait("kernel"):
                return "Jupyter Notebook"
            else:
                return "IPython Shell"

    def get_environment(self):
        return self.environment

    def get_config(self):
        config_info = {
            "Dysh Version": dysh_version,
            "Python Version": sys.version,
            "Working Directory": Path.cwd(),
            "Python Executable": Path(sys.executable),
        }
        return config_info

    def info(self):
        print("Environment: ", self.get_environment())
        print("Configuration:")
        for key, value in self.get_config().items():
            print(f"{key}: {value}")
