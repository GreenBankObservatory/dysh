import os
from pathlib import Path
from rich.align import Align
from rich.table import Table

from dysh.config.info import SystemInfo


class DyshEnvironment(SystemInfo):
    def __init__(self):
        super().__init__()
        self.dysh_base = Path(__file__).resolve().parent.parent.parent
        self.data_paths = {}
        self.get_environment()

    def get_environment(self):
        # [TODO] There's gotta be a better way to do this lol
        self.HOME = os.getenv("HOME")

        self.XDG_CONFIG_HOME = os.getenv("XDG_CONFIG_HOME")
        if self.XDG_CONFIG_HOME is None:
            self.XDG_CONFIG_HOME = os.path.join(self.HOME, ".config/")
            os.environ["XDG_CONFIG_HOME"] = self.XDG_CONFIG_HOME
            # print(f"No XDG_CONFIG_HOME found. Setting to {self.XDG_CONFIG_HOME}")

        self.DYSH_CONFIG = os.getenv("DYSH_CONFIG")
        if self.DYSH_CONFIG is None:
            self.DYSH_CONFIG = os.path.join(self.XDG_CONFIG_HOME, "dysh/")
            os.environ["DYSH_CONFIG"] = self.DYSH_CONFIG
            # print(f"No DYSH_CONFIG found. Setting to {self.DYSH_CONFIG}")

    def add_data_path(self, dir, dir_type):
        if dir not in self.data_paths.keys():
            self.data_paths[dir] = {"type": dir_type}
        else:
            print("Path already in data_dirs")

    def load_gbo_paths(self):
        self.add_data_path("/home/gbtdata", "local")

    def load_web_paths(self):
        self.add_data_path("https://www.gb.nrao.edu/dysh/example_data", "web")

    def show_data_paths(self):
        path_table = Table(title="Data Paths")
        path_table.add_column("Path", justify="left", style="green", no_wrap=True)
        path_table.add_column("Type", justify="center", style="magenta")
        for dp in self.data_paths.keys():
            path_table.add_row(dp, self.data_paths[dp]["type"])
        path_table = Align.center(path_table, vertical="middle")
        print(path_table)
