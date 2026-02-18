

*********************
dysh for GBTIDL Users
*********************


Some key differences
====================

* No longer tied to a single input and output file, can load as many SDFITS files as RAM allows

* No longer tied to a finite number of data containers, you can have as many ScanBlock/Spectrum objects as RAM allows

* Python allows for well-trackable local and global variables





Rough equivalents
=================

In the tables below, it is assumed you have executed the following commands in Python:

``from dysh.fits.gbtfitsload import GBTFITSLoad``




File I/O and Metadata Operations
--------------------------------


 .. csv-table::
    :file: files/FileIO.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10



Calibration and Data Retrieval
------------------------------


 .. csv-table::
    :file: files/Calibration.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10


Spectrum Operations
-------------------


 .. csv-table::
    :file: files/Spectrum_ops.csv
    :header-rows: 1
    :class: longtable
    :widths: 10 15 10


Side-by-side Examples
=====================

The following are examples in GBTIDL and dysh that produce equivalent results.

OTF Mapping
-----------

This example is based on the observations under TGBT17A_506_11, the same observations used for the :doc:`on-the-fly data reduction section of the users guide </users_guide/on_the_fly>`.

.. tab:: dysh

    .. code-block:: python

        sdfits = GBTFITSLoad("TGBT17A_506_11.raw.vegas")
        scan_block = sdfits.getsigref(scan=list(range(14,27)), ref=27, fdnum=0, ifnum=0, plnum=0)
        scan_block.write("dysh.fits")

.. tab:: GBTIDL

    .. code-block:: idl

        filein,"TGBT17A_506_11.raw.vegas"
        fileout,"gbtidl.fits"
        for s=14,26,1 do begin &$
            getsigref,s,27,/avgref,/keepints &$
            keep &$
        endfor

Extragalactic 21 cm Observations using Position Switching
---------------------------------------------------------

This example was borrowed from `GBTdocs position switched tutorial <https://gbtdocs.readthedocs.io/en/latest/tutorials/hi_psw_tutorial.html#data-reduction-scripted>`_.
A more thorough dysh example using the same data can be found in :doc:`example reduction of an HI survey </users_guide/hi_survey>`.
The GBTIDL version of this example can only be run with a display.

.. tab:: dysh

    .. code-block:: python

        sdfits = GBTFITSLoad("/home/astro-util/HIsurvey/Session02")

        tp0 = sdfits.gettp(scan=299, ifnum=0, plnum=0, fdnum=0).timeaverage()
        tsys0 = tp0.meta["TSYS"]
        tp1 = sdfits.gettp(scan=299, ifnum=0, plnum=1, fdnum=0).timeaverage()
        tsys1 = tp1.meta["TSYS"]

        sigref0a = sdfits.getsigref(scan=[296], ref=295, ifnum=0, fdnum=0, plnum=0, t_sys=tsys0, units="flux", zenith_opacity=0.08).timeaverage()
        sigref0b = sdfits.getsigref(scan=[298], ref=297, ifnum=0, fdnum=0, plnum=0, t_sys=tsys0, units="flux", zenith_opacity=0.08).timeaverage()
        sigref0 = sigref0a.average(sigref0b)
        sigref1a = sdfits.getsigref(scan=[296], ref=295, ifnum=0, fdnum=0, plnum=1, t_sys=tsys1, units="flux", zenith_opacity=0.08).timeaverage()
        sigref1b = sdfits.getsigref(scan=[298], ref=297, ifnum=0, fdnum=0, plnum=1, t_sys=tsys1, units="flux", zenith_opacity=0.08).timeaverage()
        sigref1 = sigref1a.average(sigref1b)
        sigref = sigref0.average(sigref1)

        sigref_smo = sigref.smooth(kernel="gauss", width=100, decimate=0)

        region = [[1.402*u.GHz, 1.4045*u.GHz],
                  [1.40506*u.GHz, 1.4054*u.GHz],
                  [1.4072*u.GHz, 1.4115*u.GHz]]

        sigref_smo.baseline(1, model="poly", include=region, remove=True)

        stats_b = sigref_smo[2000*u.km/u.s:2500*u.km/u.s].stats()
        stats_r = sigref_smo[3500*u.km/u.s:4000*u.km/u.s].stats()

        rms = (stats_b["rms"] + stats_r["rms"])/2.

        cog = sigref_smo.cog(bchan=100, echan=200)


.. tab:: GBTIDL

    .. code-block:: idl

        dirin,'/home/astro-util/HIsurvey/Session02'

        freeze      ; keep from plotting on the screen, to speed up the processing

        gettp, 299, plnum=0, /quiet
        tsys0 = !g.s.tsys
        gettp, 299, plnum=1, /quiet
        tsys1 = !g.s.tsys

        for i=295,297,2 do begin &$
          for p=0,1 do begin &$
            if (p eq 0) then tsys=tsys0[0] else tsys=tsys1[0] &$
            getsigref, i+1, i, plnum=p, tsys=tsys, unit='Jy', /quiet &$
            accum &$
          endfor &$
        endfor
        ave

        gsmooth, 100, /decimate

        setxunit, 'GHz'
        setx, 1.401, 1.412

        show

        region = [1.402, 1.4045, 1.40506, 1.4054, 1.4072, 1.4115]
        reg1 = xtochan(region)
        Nregion, reg1
        nfit, 1
        bshape
        baseline

        velo
        show
        stats, 2000, 2500, ret=mystats
        rms1 = mystats.rms
        stats, 3500, 4000, ret=mystats
        rms2 = mystats.rms

        rms = (rms1 + rms2) / 2

        gmeasure, 1, 0.5, brange=2815, erange=3142, rms=rms
        gmeasure, 1, 0.2, brange=2815, erange=3142, rms=rms
