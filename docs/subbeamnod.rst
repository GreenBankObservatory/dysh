***************
SubBeamNod Data 
***************

Calibrating SubBeamNod Data
===========================


SubBeamNod data, where a single pixel in a pixel camera is used to do a beam nod, is retrieved using :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.subbeamnod` which returns a :class:`~dysh.spectra.spectra.Spectrum` object.  First, import the relevant module::

    >>> from dysh.fits.gbtfitsload import GBTFITSLoad

..  (TODO need to replace fixed path with get_example_data() and explanation thereof)::

Then load your SDFITS file containing SubBeamNod data. In this example, we use a 
`GBT SDFITS file downloadable from GBO <http://www.gb.nrao.edu/dysh/example_data/subbeamnod-Ka/data/TRCO_230413_Ka.raw.vegas/TRCO_230413_Ka.raw.vegas.A.fits'>_``

    >>> f = 'TRCO_230413_Ka.raw.vegas.A.fits'
    >>> sdfits = GBTFITSLoad(f)

The returned `sdfits` can be probed for information::

    >>> sdfits.info()
    Filename: /data/gbt/examples/subbeamnod-Ka/data/TRCO_230413_Ka.raw.vegas/TRCO_230413_Ka.raw.vegas.A.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU      12   ()      
      1  SINGLE DISH    1 BinTableHDU    245   5280R x 74C   ['32A', '1D', '22A', '1D', '1D', '1D', '1024E', '16A', '6A', '8A', '1D', '1D', '1D', '4A', '1D', '4A', '1D', '1I', '32A', '32A', '1J', '32A', '16A', '1E', '8A', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '8A', '1D', '1D', '12A', '1I', '1I', '1D', '1D', '1I', '1A', '1I', '1I', '16A', '16A', '1J', '1J', '22A', '1D', '1D', '1I', '1A', '1D', '1E', '1D', '1D', '1D', '1D', '1D', '1A', '1A', '8A', '1E', '1E', '16A', '1I', '1I', '1I']   

You can also print a concise (or verbose if you choose `verbose=True`) summary :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.summary` of the data::

    >>> sdfits.summary()
            SCAN     OBJECT VELOCITY        PROC PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
        0   32.0  1256-0547      0.0         Nod      1.0     26.5    26.5    1     2    30      2  160.975324  43.884984
        1   33.0  1256-0547      0.0         Nod      2.0     26.5    26.5    1     2    30      2  161.174093  43.928449
        2   34.0  1256-0547      0.0         Nod      1.0     30.5    30.5    1     2    30      2  161.589629  44.000491
        3   35.0  1256-0547      0.0         Nod      2.0     30.5    30.5    1     2    30      2  161.783395  44.041622
        4   36.0  1256-0547      0.0     Unknown      0.0     0.75    0.75    1     2   120      2  162.124052  44.100404
        5   37.0  1256-0547      0.0         Nod      1.0     34.5    34.5    1     2    30      2  162.611075  44.183661
        6   38.0  1256-0547      0.0         Nod      2.0     34.5    34.5    1     2    30      2  162.896506  44.237997
        7   39.0  1256-0547      0.0         Nod      1.0     37.5    37.5    1     2    30      2  163.333508  44.306385
        8   40.0  1256-0547      0.0         Nod      2.0     37.5    37.5    1     2    30      2  163.529285  44.343704
        9   41.0  1256-0547      0.0         Nod      1.0     30.5    30.5    1     2    30      2  164.941425  44.559629
        10  42.0  1256-0547      0.0         Nod      2.0     30.5    30.5    1     2    30      2  165.139436  44.593378
        11  43.0  1256-0547      0.0  SubBeamNod      1.0     30.5    30.5    1     2   120      2  165.469522  44.639023
        12  44.0  1256-0547      0.0         Nod      1.0     30.5    30.5    1     2    30      2   166.48287  44.776997
        13  45.0  1256-0547      0.0         Nod      2.0     30.5    30.5    1     2    30      2  166.688378  44.808119
        14  46.0  1256-0547      0.0  SubBeamNod      1.0     30.5    30.5    1     2   120      2  167.026583  44.849753
        15  52.0  1256-0547      0.0         Nod      1.0     30.5    30.5    1     2    30      2  169.972904  45.179358
        16  53.0  1256-0547      0.0         Nod      2.0     30.5    30.5    1     2    30      2  170.175815  45.201877
        17  54.0  1256-0547      0.0  SubBeamNod      1.0     30.5    30.5    1     2   120      2  170.518885  45.232575

The SubBeamNod scans are 43, 46, and 54.  Retrieve and calibrate a SubBeamNod scan, then plot it::

    >>> sbn = sdfits.subbeamnod(scan=43,fdnum=1,ifnum=0,weights='tsys')
    >>> sbn.plot()

