*********************************
dysh Interactive Plotter Tutorial
*********************************

In this tutorial, we will walk through the features of the interactive plotter offered by dysh.
There are two types of interactive plots: waterfall/spectrogram plots, and 1-dimensional spectrum plots. 



=========================
ScanBlock Waterfall plots
=========================

In this section, we describe how to create and interact with waterfall plots.
First, we download data from an HI survey and open it with dysh.

.. code-block::

    from dysh.fits.gbtfitsload import GBTFITSLoad
    from pathlib import Path
    from dysh.util.download import from_url


.. code-block::

    url = "http://www.gb.nrao.edu/dysh/example_data/rfi-L/data/AGBT17A_404_01.raw.vegas/AGBT17A_404_01.raw.vegas.A.fits"
    savepath = Path.cwd() / "data"
    filename = from_url(url, savepath)
    sdfits = GBTFITSLoad(filename)


Now, we grab a position-switched scan with GPS interference but we don't time average it, leaving it as a full ScanBlock with one PSScan.
Waterfall plots can be created using single Scans or Scanblocks that contain many scans. The process is identical.


.. code-block::

    ps = sdfits.getps(scan=19, plnum=0, fdnum=0, ifnum=0)
    psplt = ps.plot()

.. image:: files/sb_plot1.png

Since we keep a reference to the plotter object `psplt`, we can do several things at the command line to change its aspects.
Let's change the lower limit of the color scale to see the clean data and the extent of the GPS RFI a little better.

.. code-block::

    psplt.set_clim(vmin = -2)

.. image:: files/sb_plot2.png

Other CLI options include `set_interpolation()`, `set_cmap()`, and `set_norm()`.
All of the kwargs associated with these functions can be added to the initial plot instantiation as well.
See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html>`_ for more details on these kwargs.

You can also do waterfall plots of scanblocks with more than one scan.

.. code-block::

    tp = sdfits.gettp(scan=[19,20], plnum=0, fdnum=0, ifnum=0)
    tpplt = tp.plot(interpolation='gaussian',cmap='hot',norm='log',vmin=3e8,vmax=4e8)

Here we see clearly that the GPS RFI turns on towards the end of the second scan.
The X-axis denotes the scan number on the bottom, and the integration number along the tick marks.
The integration number resets to 0 at the beginning of each scan.

.. image:: files/sb_plot3.png


=========================
Spectrum plots
=========================


.. code-block::

    from dysh.fits.gbtfitsload import GBTFITSLoad
    from pathlib import Path
    from dysh.util.download import from_url


.. code-block::

    url = "http://www.gb.nrao.edu/dysh/example_data/rfi-L/data/AGBT17A_404_01.raw.vegas/AGBT17A_404_01.raw.vegas.A.fits"
    savepath = Path.cwd() / "data"
    filename = from_url(url, savepath)
    sdfits = GBTFITSLoad(filename)


Now grab a position-switched scan with GPS interference, average it, and smooth
with a 7-channel wide boxcar kernel. Start the interactive plotter with the `plot()` command.
`Spectrum.plot`` supports a subset of the arguments offered by `matplotlib.pyplot.plot`. These are: `xmin`, `xmax`,
`ymin`, `ymax`, `xlabel`, `ylabel`, `xaxis_unit`, `yaxis_unit`, `grid`, `linewidth`, `color`, and `title`.
If you do not wish to have the interactive plotter and prefer a static plot, add `interactive=False`
to the arguments.


.. code-block::

    ps = sdfits.getps(scan=19, plnum=0, fdnum=0, ifnum=0).timeaverage()
    ps.plot()


.. image:: files/iplotter-2.png

Use the matplotlib Zoom button (red) to Zoom into the baseline.


.. image:: files/iplotter-3.png

We will use a first order baseline. Uncheck the Zoom button (red) to re-enable region selection,
and go back to the original zoom level with the Home button (green).


.. image:: files/iplotter4.gif

Create a region over the clean part of the spectrum by holding down left-click with your mouse
and dragging across the screen.

.. image:: files/iplotter5.png

You can clear all regions with the "Clear Regions" (red) button.

.. image:: files/iplotter6.gif

You can click on a region to select it, allowing you to drag it to a different place in the spectrum or
change its size. Currently, regions can only be resized by clicking near the edge from within the region.
Regions are allowed to overlap. Once a region is selected,
you can delete it with the "Delete Region" button.

.. admonition:: Note

    All buttons on the plot will disappear when using the `ps.savefig()` command. However, they will still
    remain in screenshots.


.. image:: files/iplotter7.png

With one large region selected, we can try a baseline. Use the `remove=False` option to
just plot the baseline fit without removing it.

.. code-block::

    ps.baseline(1,include=ps.get_selected_regions(),remove=False)


The fit looks good. Let's subtract it out and save the figure. See that the buttons
have disappeared. The white space at the top is reserved for buttons that will be implemented
in the future and will be addressed. We can check the stats before and after to see they've improved.

.. code-block::

    ps.get_selected_regions()
    [(210,2871)]

    print(f"{ps[210:2871].stats()['mean']:.4f}, {ps[210:2871].stats()['rms']:.4f}")
    -0.0780 K, 0.0204 K

    ps.baseline(1,include=ps.get_selected_regions(),remove=True)

    print(f"{ps[210:2871].stats()['mean']:.4f}, {ps[210:2871].stats()['rms']:.4f}")
    -0.0000 K, 0.0204 K

    ps.savefig('my_spectrum.png')

.. image:: files/iplotter8.png
