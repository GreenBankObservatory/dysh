import sys, time
from IPython import display
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QApplication, QMainWindow
from plotly.graph_objects import Figure, Scatter
import plotly
import numpy as np
from dash import Dash, html, dash_table, dcc, callback, Output, Input
import dash_bootstrap_components as dbc

from dysh.plot.renderer import Renderer

class SpecPlot:

    def __init__(self, spectrum, **kwargs):
        # self.reset()
        self._spectrum = spectrum
        # self._set_xaxis_info()
        # self._plot_kwargs.update(kwargs)
        # self._plt = plt
        # self._figure = None
        # self._axis = None
        # self._title = self._plot_kwargs["title"]

        self.detect_renderer()
        self.get_x()
        self.get_y()
        self.get_mask()
        self.make_fig()
        # self.add_dropdowns()
        self.plot_spectrum()
        self._init_dash()
        self.build_html()

    def _init_dash(self):
        self.app = Dash()
        self.app.layout = html.Div([
            dcc.Textarea(
                id='title',
                value='Dysh Plotter',
                ),
            dbc.Row([
                dcc.Graph(figure={}, id='graph-placeholder'),
                dcc.Checklist(
                    id='plot-checklist',
                    options=[
                        {'label': 'Gridlines', 'value': 'gridlines'},
                    ],
                    value=['gridlines']
                )
            ]),
        ])

    def detect_renderer(self):
        self.renderer = Renderer()
        self.renderer.info()

    def get_x(self):
        self.x = self._spectrum.spectral_axis
    
    def get_y(self):
        self.y = self._spectrum.flux

    def get_mask(self):
        self.mask = self._spectrum.mask

    def make_toggle_gridlines(self):
        self.toggle_grid = dcc.Checklist(
            options=[
                {'label': 'Grid', 'value': 'grid'},
                ],
            value=['Montreal']
            )
        
    def make_fig(self):
        self.fig = Figure()
        self.update_fig()

    # @callback(
    #     Output(component_id='graph-placeholder', component_property='figure'),
    #     Input(component_id='plot-checklist', component_property='value')
    #     )
    def update_fig(self):
        self.set_title()
        self.set_xlabel()
        self.set_ylabel()
        self.update_xaxis()
        self.update_yaxis()
        self.update_background()
        self.plot_unmasked()
        self.plot_masked()    

    def set_title(self):
        self.fig.update_layout(
            title=dict(
                text="Sample Spectrum",
                x=0.5,
                ),
        )

    def set_xlabel(self):
        self.fig.update_layout(
            xaxis_title=dict(text="Frequency")
        )

    def set_ylabel(self):
        self.fig.update_layout(
            yaxis_title=dict(text="Flux")
        )

    def update_xaxis(self):
        self.fig.update_xaxes(
            ticks='outside',
            showline=True,
            linecolor='white',
            gridcolor='lightgrey'
            )
        
    def update_yaxis(self):
        self.fig.update_yaxes(
            ticks='outside',
            showline=True,
            linecolor='white',
            gridcolor='lightgrey'
            )

    def update_background(self):
        self.fig.update_layout(
            paper_bgcolor='black',
            plot_bgcolor='black',
            font=dict(color='white'),
            )

    def plot_unmasked(self):
        self.trace_unmasked = dict(
            name="unmasked",
            type="scatter",
            #color='green',
            x=self.x[~self.mask],
            y=self.y[~self.mask],
            mode='lines',
            marker=dict(color='green',
                        #size=10, 
                        showscale=False)
                        )
        self.fig.add_trace(self.trace_unmasked)

    def plot_masked(self):
        self.trace_masked = dict(
            name="masked",
            type="scatter",
            x=self.x[self.mask],
            y=self.y[self.mask],
            mode='lines',
            marker=dict(color='red',
                        #size=10, 
                        showscale=False)
                        )    
        self.fig.add_trace(self.trace_masked)

    def plot_spectrum(self):
        self.plot_html = plotly.offline.plot(self.fig, output_type='div', include_plotlyjs='cdn') 

    def show(self):
        pass
        # self.app.run(debug=True)
    
    def build_html(self):
        self.html = '<html><body>'
        self.html += self.plot_html
        self.html += '</body></html>'

    def add_dropdowns(self):
        self.fig.update_layout(
            updatemenus=[
                dict(
                    buttons=list([
                        dict(
                            args=["colorscale", "Viridis"],
                            label="Viridis",
                            method="restyle"
                        ),
                        dict(
                            args=["colorscale", "Cividis"],
                            label="Cividis",
                            method="restyle"
                        ),
                        dict(
                            args=["colorscale", "Blues"],
                            label="Blues",
                            method="restyle"
                        ),
                        dict(
                            args=["colorscale", "Greens"],
                            label="Greens",
                            method="restyle"
                        ),
                    ]),
                )
            ]
        )

        self.fig.update_layout(
            annotations=[
                dict(text="Frame", x=0, xref="paper", y=1.06, yref="paper",
                                    align="left", showarrow=False),
            ])

class SpecPlotWindow(QMainWindow):

    def __init__(self, spectrum, **kwargs):

        super(SpecPlotWindow, self).__init__()

        # some example data
        my_plot = SpecPlot(spectrum, **kwargs)
        my_plot.show()
        #breakpoint()

        # we create an instance of QWebEngineView and set the html code
        plot_widget = QWebEngineView()
        plot_widget.setHtml(my_plot.html)

        # set the QWebEngineView instance as main widget
        self.setCentralWidget(plot_widget)

def make_window():
    app = QApplication([])
    window = SpecPlotWindow(ta)
    window.show()
    app.exec_()
    # while True:
    #     try:
    #         print("it's open")
    #         display.display(window)
    #         display.clear_output(wait=True)
    #         time.sleep(1)
    #     except KeyboardInterrupt:
    #         break

if __name__ == '__main__':
    from dysh.fits.gbtfitsload import GBTFITSLoad
    filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
    sdfits = GBTFITSLoad(filename)
    psscan = sdfits.getps(scan=152, ifnum=0, plnum=0)
    ta = psscan.timeaverage(weights='tsys')
    ta.mask[0:300] = True
    # breakpoint()

    make_window()
