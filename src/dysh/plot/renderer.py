import sys


class Renderer:
    def __init__(self):
        self.render_type = None
        self.detect()

    def detect(self):
        self.detect_notebook()
        self.detect_interactive()

    def info(self):
        print(f"In notebook: {self.in_notebook}")
        print(f"Interactive: {self.in_interactive}")

    def detect_notebook(self):
        """
        Returns ``True`` if the module is running in IPython kernel,
        ``False`` if in IPython shell or other Python shell.
        """
        self.in_notebook = "ipykernel" in sys.modules
        if self.in_notebook:
            self.render_type = "notebook"

    def detect_interactive(self):
        """
        Returns ``True`` if the module is running in an interactive setting,
        ``False`` otherwise.
        """
        self.in_interactive = bool(getattr(sys, "ps1", sys.flags.interactive))
        if self.in_interactive:
            self.render_type = "iPython"
