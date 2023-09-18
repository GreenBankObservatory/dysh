from bs4 import BeautifulSoup
import os, wget, socket, requests, subprocess
from getpass import getuser
from pathlib import Path
from rich.text import Text
from rich.tree import Tree

from dysh.config.environment import DyshEnvironment
from dysh.config.rich_theme import DyshRichConsole
from dysh.util.messages import SystemMessages

dysh_env = DyshEnvironment()
dysh_env.load_gbo_paths()
dysh_env.load_web_paths()


class DirectoryModel:
    def __init__(self, model_type, base_dir):
        self.model_type = model_type
        self.base_dir = base_dir
        self.console = DyshRichConsole()
        self.tree = None

    def print_tree(self):
        if self.tree is None:
            self.get_tree()
        self.console.print(self.tree)

    def get_tree(self, max_depth=10):
        self.max_depth = max_depth
        if self.model_type == "web":
            self.get_tree_web()

    def get_tree_web(self):
        # [TODO] This can print a LOT. Maybe summarize more?
        self.tree = Tree(
            f":open_file_folder: [link {self.base_dir}]{self.base_dir}",
            guide_style="bold bright_blue",
        )
        self._get_tree_web(self.tree, self.base_dir)

    def _get_tree_web(self, tree, url, current_depth=0):
        ignore_starts = ["?", "/"]
        if current_depth <= self.max_depth:
            current_depth += 1
            page = requests.get(url).text
            soup = BeautifulSoup(page, "html.parser")
            links = [node.get("href") for node in soup.find_all("a")]
            links.sort()
            for li in links:
                continue_tree = True
                for istart in ignore_starts:
                    if li.startswith(istart):
                        continue_tree = False
                if continue_tree:
                    full_url = os.path.join(url, li)
                    if li.endswith("/"):
                        branch = tree.add(f"[bold magenta]:open_file_folder: [link {full_url}]{li}")
                        self._get_tree_web(branch, full_url, current_depth)
                    else:
                        text_filename = Text(li, "green")
                        text_filename.stylize(f"link {full_url}")
                        icon = "🐍 " if li.endswith(".py") else "📄 "
                        tree.add(Text(icon) + text_filename)


url = "https://www.gb.nrao.edu/dysh/example_data/"
web_examples = DirectoryModel("web", url)
web_examples.print_tree()