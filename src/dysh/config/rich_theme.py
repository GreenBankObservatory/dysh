from rich.console import Console
from rich.theme import Theme
from rich.highlighter import RegexHighlighter


class DyshRichHighlighter(RegexHighlighter):
    """Apply style to certain text"""

    # [TODO] this doesn't work yet. why must regex hurt me in this way

    base_style = "file_extension."
    highlights = [
        r"(?P<fits>\.fits)",
        r"(?P<flag>\.flag)",
        r"(?P<fits>\.flag)",
    ]


class DyshRichConsole(Console):
    def __init__(self):
        self._init_theme()
        self._init_highlighter()
        Console.__init__(self, highlighter=self.highlighter, theme=self.custom_theme)

    def _init_theme(self):
        self.custom_theme = Theme(
            {
                "file_extension.fits": "magenta",
                "file_extension.flag": "bold red",
                "file_extension.index": "dim cyan",
            }
        )

    def _init_highlighter(self):
        self.highlighter = DyshRichHighlighter()
