import sys, os
import json, pathlib
from rich.align import Align
from rich.filesize import decimal
from rich.markup import escape
from rich.text import Text
from rich.tree import Tree
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint
from dysh.config.rich_theme import DyshRichConsole

sys.path.insert(0, os.path.abspath("."))


class SystemMessages:
    """Messages about the system"""

    def __init__(self):
        pass

    def directory_tree(self, dir, max_depth=4):
        tree = Tree(
            f":open_file_folder: [link file://{dir}]{dir}",
            guide_style="bold bright_blue",
        )
        self._directory_tree(pathlib.Path(dir), tree, max_depth)
        DyshRichConsole.print(tree)

    def _directory_tree(self, directory, tree, max_depth, current_depth=0):
        """Recursively build a Tree with directory contents"""

        if current_depth <= max_depth:
            current_depth += 1
            # Sort dirs first,  then sort by filename
            paths = sorted(
                pathlib.Path(directory).iterdir(),
                key=lambda path: (path.is_file(), path.name.lower()),
            )

            for path in paths:
                # Ignore certain files
                ignore_starts = [".", "_"]
                continue_tree = True
                for istart in ignore_starts:
                    if path.name.startswith(istart):
                        continue_tree = False

                if continue_tree:
                    if path.is_dir():
                        style = ""
                        branch = tree.add(
                            f"[bold magenta]:open_file_folder: [link file://{path}]{escape(path.name)}",
                            style=style,
                            guide_style=style,
                        )
                        self._directory_tree(path, branch, max_depth, current_depth)
                    else:
                        text_filename = Text(path.name, "green")
                        text_filename.highlight_regex(r"\..*$", "bold red")
                        text_filename.stylize(f"link file://{path}")
                        file_size = path.stat().st_size
                        text_filename.append(f" ({decimal(file_size)})", "blue")
                        icon = "🐍 " if path.suffix == ".py" else "📄 "
                        tree.add(Text(icon) + text_filename)


class GBTInfoMessages:
    """Messages about the GBT"""

    def __init__(self):
        pass

    def RXTable(self):
        f = open("/Users/victoriacatlett/Documents/Code/repos/github/dysh/gui/static/gbt/gbt-rx.json")
        rx_dict = json.load(f)

        rx_table = Table(title="GBT Receivers")

        rx_table.add_column("Name", justify="right", style="cyan", no_wrap=True)
        rx_table.add_column("Frequencies", style="magenta")

        for rx in rx_dict:
            rx_name = rx.split(" ")[0]
            rx_freqs = rx.split("(")[1][0:-1]
            rx_table.add_row(rx_name, rx_freqs)

        rx_table = Align.center(rx_table, vertical="middle")
        # DyshRichConsole.print(rx_table)
        rprint(rx_table)


class FriendlyMessages:
    """General friendly stuff to print"""

    def hello():
        w_message = "Welcome to Dysh"
        w_panel = Panel(w_message, padding=2)
        w_panel = Align.center(w_panel, vertical="middle")
        # DyshRichConsole.print(w_panel)
        rprint(w_panel)

    def goodbye():
        bye_message = "Goodbye!"
        bye_panel = Panel(bye_message, padding=2)
        bye_panel = Align.center(bye_panel, vertical="middle")
        # DyshRichConsole.print(bye_panel)
        rprint(bye_panel)