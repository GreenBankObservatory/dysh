"""
Base class for dysh GUIs
"""


class BaseGUI:
    def connect_buttons(self, plotbase):
        """
        Connect the GUI button(s) callbacks.

        Parameters
        ----------
        plotbase : `~dysh.plot.plotbase.PlotBase`
            Plotting class that defines the functions to be used in the callbacks.
        """
        return

    def disconnect_buttons(self):
        """
        Disconnect the GUI button(s) callbacks.
        """
        return

    def is_window_alive(self):
        """
        Is the GUI window alive.

        Returns
        -------
            True if the GUI window is still alive.
        """
        return True

    def show(self):
        """
        Show the GUI contents.
        """
        return
