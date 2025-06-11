*********************************
dysh Interactive Plotter Tutorial
*********************************




In this tutorial, we will walk through the features of the interactive plotter offered by dysh.
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


Now grab a position-switched scan with GPS interference, average it, and smooth
with a 7-channel wide boxcar kernel. Start the interactive plotter with the `plot()` command. This command 
supports all kwargs offered by the standard `matplotlib.pyplot.plot()` command.
If you do not wish to have the interactive plotter and prefer a static plot, add `interactive=False` to the arguments.


.. code-block::

    ps_sm = sdfits.getps(scan=19, plnum=0, fdnum=0, ifnum=0).timeaverage().smooth('boxcar',7)
    ps_sm.plot()








