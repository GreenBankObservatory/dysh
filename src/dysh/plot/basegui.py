"""
Base class for dysh GUIs
"""


class BaseGUI:
    def connect_buttons(self, plotbase):
        return

    def disconnect_buttons(self):
        return

    def is_window_alive(self):
        return True

    def show(self):
        return
