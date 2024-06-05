#%%
from dysh.fits.gbtfitsload import GBTFITSLoad
from bokeh.io import show, reset_output, output_notebook
from bokeh.layouts import column
from bokeh.models import ColorPicker
from bokeh.plotting import figure

reset_output()
output_notebook()
#%%
filename = "E:/Code/GitHub/Work/forks/dysh/test_file.fits"
sdfits = GBTFITSLoad(filename)
psscan = sdfits.getps(152, ifnum=0, plnum=0)
ta = psscan.timeaverage(weights='tsys')
#%%
print(dir(ta))
# %%
plot = figure(x_range=(0, 1), y_range=(0, 1), width=350, height=350)
line = plot.line(x=(0,1), y=(0,1), color="black", line_width=4)

picker = ColorPicker(title="Line Color")
picker.js_link('color', line.glyph, 'line_color')

show(column(plot, picker))