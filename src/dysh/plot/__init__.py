"""Classes and functions for plotting spectra and SDFITS data"""

# Interactive plotting.
import matplotlib.pyplot as plt

plt.ion()

__all__ = ["specplot"]
from .core import *
from .specplot import *
from .scanplot import *
from .vegasplot import *
