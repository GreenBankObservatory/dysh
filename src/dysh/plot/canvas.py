import matplotlib.pyplot as plt
import numpy as np
import ipywidgets as widgets
from IPython.display import display

# Generate sample data
x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)
print(x)
print(y)
# Create a Matplotlib figure and a set of axes
fig, ax = plt.subplots(1)
# Plot the data
line, = ax.plot(x, y)
# Define a function to update the plot based on the selected value
def update_plot(value):
    # Update the y-values based on the selected value
    if value == 'sin(x)':
        line.set_ydata(np.sin(x))
    elif value == 'cos(x)':
        line.set_ydata(np.cos(x))
    else:
        line.set_ydata(np.tan(x))
    # Redraw the plot
    fig.canvas.draw()
# Create a drop-down widget
dropdown = widgets.Dropdown(
    options=['sin(x)', 'cos(x)', 'tan(x)'],
    value='sin(x)',
    description='Function:'
    )
# Register the update_plot function as the callback for the dropdown widget
dropdown.observe(update_plot, 'value')
# Display the dropdown and the plot
display(dropdown)
plt.show()