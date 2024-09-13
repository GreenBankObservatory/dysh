from dysh.plot.canvas import SpectrumPlot
from dysh.plot.renderer import EnvironmentInfo


def plot_spectrum(spectrum=None, display=True, interactive=False):
    env_info = EnvironmentInfo()
    renderer = env_info.get_environment()

    if renderer == "Python Script":
        return plot_spectrum_script(spectrum, display, interactive)
    elif renderer == "Jupyter Notebook":
        return plot_spectrum_jupyter(spectrum, display, interactive)
    elif renderer == "IPython Shell":
        return plot_spectrum_ipython(spectrum, display, interactive)
    else:
        print("Unknown Python configuration")
        return None


def plot_spectrum_script(spectrum, display=True, interactive=False):
    my_plot = SpectrumPlot(spectrum, interactive)
    if display:
        my_plot.show()
    return my_plot


def plot_spectrum_jupyter(spectrum, display=True, interactive=False):
    my_plot = SpectrumPlot(spectrum, interactive)
    return my_plot


def plot_spectrum_ipython(spectrum, display=True, interactive=False):
    my_plot = SpectrumPlot(spectrum, interactive)
    if display:
        my_plot.show()
    return my_plot
