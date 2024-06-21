from IPython import get_ipython


class Renderer:
    def __init__(self):
        self.render_type = None
        self.detect()

    def detect(self):
        """Detect the renderer needed"""
        ip = get_ipython()

        if ip is None:
            self.render_type = "script"
        else:
            if ip.has_trait("kernel"):
                self.render_type = "notebook"
            else:
                self.render_type = "ipython"

    def info(self):
        """Print out the detected renderer"""
        print(f"Renderer: {self.render_type}")
