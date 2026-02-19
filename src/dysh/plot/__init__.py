"""Classes and functions for plotting spectra and SDFITS data"""

__all__ = ["core", "specplot", "scanplot", "vegasplot"]


def switch_frontend(frontend):
    from dysh.plot import plotbase

    plotbase.GUI = frontend
