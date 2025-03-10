"""
Core functions for FITS/SDFITS
"""

import matplotlib.pyplot as plt
import numpy as np


def default_sdfits_columns():
    """The default column names for GBT SDFITS.

    Returns
    -------
        colnames : list
        A list of the GBT SDFITS column names
    """
    colnames = [
        "OBJECT",
        "BANDWID",
        "DATE-OBS",
        "DURATION",
        "EXPOSURE",
        "TSYS",
        "TDIM7",
        "TUNIT7",
        "CTYPE1",
        "CRVAL1",
        "CRPIX1",
        "CDELT1",
        "CTYPE2",
        "CRVAL2",
        "CTYPE3",
        "CRVAL3",
        "CRVAL4",
        "OBSERVER",
        "OBSID",
        "SCAN",
        "OBSMODE",
        "FRONTEND",
        "TCAL",
        "VELDEF",
        "VFRAME",
        "RVSYS",
        "OBSFREQ",
        "LST",
        "AZIMUTH",
        "ELEVATIO",
        "TAMBIENT",
        "PRESSURE",
        "HUMIDITY",
        "RESTFREQ",
        "DOPFREQ",
        "FREQRES",
        "EQUINOX",
        "RADESYS",
        "TRGTLONG",
        "TRGTLAT",
        "SAMPLER",
        "FEED",
        "SRFEED",
        "FEEDXOFF",
        "FEEDEOFF",
        "SUBREF_STATE",
        "SIDEBAND",
        "PROCSEQN",
        "PROCSIZE",
        "PROCSCAN",
        "PROCTYPE",
        "LASTON",
        "LASTOFF",
        "TIMESTAMP",
        "QD_XEL",
        "QD_EL",
        "QD_BAD",
        "QD_METHOD",
        "VELOCITY",
        "FOFFREF1",
        "ZEROCHAN",
        "ADCSAMPF",
        "VSPDELT",
        "VSPRVAL",
        "VSPRPIX",
        "SIG",
        "CAL",
        "CALTYPE",
        "TWARM",
        "TCOLD",
        "CALPOSITION",
        "BACKEND",
        "PROJID",
        "TELESCOP",
        "SITELONG",
        "SITELAT",
        "SITEELEV",
        "IFNUM",
        "PLNUM",
        "FDNUM",
        "INT",
        "INTNUM",  # not all SDFITS files have INT, so we always create INTNUM
        "NSAVE",
        # The following are added by the GBTFITSLoad constructor.
        # Arguable whether they should be included or not.
        "HDU",
        "BINTABLE",
        "ROW",
        "PROC",
        "OBSTYPE",
        "SUBOBSMODE",
    ]
    return colnames


def mean_data(data, fedge=0.1, nedge=None, median=False):
    """
    special mean of data to exclude the edges like mean_tsys(), with
    an option to use the median instead of the mean.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data.
    fedge : float, optional
        Fraction of edge channels to exclude at each end, a number between 0 and 1.
        If `nedge` is used, this parameter is not used.
        Default: 0.1, meaning the central 80% bandwidth is used
    nedge : int, optional
        Number of edge channels to exclude. nedge cannot be 0.
        Default: None, meaning use `fedge`
    median : boolean, optional
        Use the median instead of the mean.
        The default is False.

    Returns
    -------
    meandata : float

    """

    nchan = len(data)
    if nedge is None:
        nedge = int(nchan * fedge)
    chrng = slice(nedge, -(nedge - 1), 1)
    if median:
        meandata = np.nanmedian(data[chrng])
    else:
        meandata = np.nanmean(data[chrng])
    return meandata


def getbeam(sdf, debug=False):
    """
    find the two nodding beams based on on a given FDNUM, FEED
    needs PROCSCAN='BEAM1' or 'BEAM2'

    Parameters
    ----------
    sdf : `GBTFITSLoad`
        data handle, containing one or more SDFITS files specific to GBT
    debug : boolean, optional
        Add more debugging output. @todo use logger
        The default is False.

    Returns
    -------
    beams : list of two ints representing the nodding beams (0 = first beam)

    """
    # list of columns needed to differentiate and find the nodding beams
    kb = ["FEEDXOFF", "FEEDEOFF", "PROCSCAN", "FDNUM", "FEED"]
    a0 = sdf._index[kb]
    b1 = a0.loc[a0["FEEDXOFF"] == 0.0]
    b2 = b1.loc[b1["FEEDEOFF"] == 0.0]
    d1 = b2.loc[b2["PROCSCAN"] == "BEAM1"]
    d2 = b2.loc[b2["PROCSCAN"] == "BEAM2"]
    #
    if len(d1["FDNUM"].unique()) == 1 and len(d2["FDNUM"].unique()) == 1:
        beam1 = d1["FDNUM"].unique()[0]
        beam2 = d2["FDNUM"].unique()[0]
        fdnum1 = d1["FEED"].unique()[0]
        fdnum2 = d2["FEED"].unique()[0]
        if debug:
            print("beams: ", beam1, beam2, fdnum1, fdnum2)
        return [beam1, beam2]
    else:
        # try one other thing
        if len(b2["FEED"].unique()) == 2:
            print("getbeam rescued")
            b = b2["FEED"].unique() - 1
            return list(b)
        print("too many in beam1:", d1["FDNUM"].unique())
        print("too many in beam2:", d2["FDNUM"].unique())
        return []


def calseq(sdf, scan, tcold=54, fdnum=0, ifnum=0, plnum=0, freq=None, verbose=False):
    """
    This routine returns the Tsys and gain for the selected W-band channel.

    W-band receivers use a CALSEQ where during a scan three different
    observations are made: sky, cold1 and cold2, from which the
    system temperature is derived.


    Parameters
    ----------
    sdf : `GBTFITSLoad`
        data handle, containing one or more SDFITS files specific to GBT
    scan : int or list of int
        Scan number(s) where CALSEQ is expected. See sdf.summary() to find the scan number(s).
        If multiple scans are used, an average Tsys is computed.
    tcold : float, optional
        Set the cold temperature. See also freq= for an alternative computation.
        The default is 54.
    fdnum : int, optional
        Feed to be used, 0 being the first.
        The default is 0.
    ifnum : int, optional
        IF to be used, 0 being the first.
        The default is 0.
    plnum : int, optional
        Polarization to be used, 0 being the first.
        The default is 0.
    freq : float, optional
        Observing frequency if Tcold to be set different from the default:
        Tcold = 54 - 0.6*(FREQ-77)      FREQ in GHz
        The default is None.
    verbose : boolean, optional
        Add more information mimicking the GBTIDL outout of VANECAL.
        The default is False

    Returns
    -------
    tsys : float
        DESCRIPTION.
    g : float
        DESCRIPTION.

    """
    if freq is not None:
        # see eq.(13) in GBT memo 302
        tcold = 54 - 0.6 * (freq - 77)
        print(f"Warning: calseq using freq={freq} GHz and setting tcold={tcold} K")

    twarm = sdf._index["TAMBIENT"].mean()

    tp_args = {"scan": scan, "ifnum": ifnum, "plnum": plnum, "fdnum": fdnum, "calibrate": True, "cal": False}
    vsky = sdf.gettp(CALPOSITION="Observing", **tp_args).timeaverage()
    vcold1 = sdf.gettp(CALPOSITION="Cold1", **tp_args).timeaverage()
    vcold2 = sdf.gettp(CALPOSITION="Cold2", **tp_args).timeaverage()

    if fdnum == 0:
        g = (twarm - tcold) / mean_data(vcold2.data - vcold1.data)
    elif fdnum == 1:
        g = (twarm - tcold) / mean_data(vcold1.data - vcold2.data)
    else:
        print(f"Illegal fdnum={fdnum} for a CALSEQ")
        return None
    tsys = mean_data(g * vsky.data)

    if verbose:
        print(f"Twarm={twarm} Tcold={tcold}")
        print(f"IFNUM {ifnum} PLNUM {plnum} FDNUM {fdnum}")
        print(f"Tsys = {tsys}")
        print(f"Gain [K/counts] = {g}")

    return tsys, g


def vanecal(sdf, vane_sky, feeds=range(16), mode=2, tcal=None, verbose=False):
    """
    Return Tsys calibration values for all or selected beams of the ARGUS
    VANE/SKY calibration cycle.


    Parameters
    ----------
    sdf : `GBTFITSLoad`
        data handle, containing one or more SDFITS files specific to GBT
    vane_sky : list of two ints
        The first designates the VANE scan, the second the SKY scan.
        Normally the SKY scan is directly followed by the VANE scan.
        @todo if one scan given, assume sky is vane+1
    feeds : list of ints, optional
        The default is range(16), i.e. using all ARGUS beams.
    mode : int, optional
        Mode of computing. See also `mean_tsys()`
        mode=0  Do the mean before the division
        mode=1  Do the mean after the division
        mode=2  Take a median of the inverse division
        The default is 2.
    tcal : float, optional
        Tcal value for normalization. Normally obtained from the
        environment, but offsite cannot be done.
        @todo fix this, but right now it is adviced to manually enter tcal
        The default is None.
    verbose : boolean, optional
        Add more information mimicking the GBTIDL outout of VANECAL.
        The default is False

    Returns
    -------
    tsys : list of floats
        Values of Tsys for each of the `feeds` given.

    """
    vane = vane_sky[0]
    sky = vane_sky[1]
    if len(feeds) == 0:
        print("Warning, no feeds= given")
        return None
    tsys = np.zeros(len(feeds), dtype=float)
    sdf1 = sdf

    #  for VANE/CAL data usually tcal=1
    if tcal is None:
        tcal = sdf._index["TCAL"].mean()
        if tcal == 1.0:
            # until we figure this out via getatmos
            tcal = 100.0  # set to 100K for now

    i = 0
    for f in feeds:
        v = sdf1.gettp(scan=vane, fdnum=f, calibrate=True, cal=False).timeaverage()
        s = sdf1.gettp(scan=sky, fdnum=f, calibrate=True, cal=False).timeaverage()
        if mode == 0:
            mean_off = mean_data(s.data)
            mean_dif = mean_data(v.data - s.data)
            tsys[i] = tcal * mean_off / mean_dif
        elif mode == 1:
            tsys[i] = tcal / mean_data((v.data - s.data) / s.data)
        elif mode == 2:
            tsys[i] = tcal / np.nanmedian((v.data - s.data) / s.data)
        #  vanecal.pro seems to do    tcal / median( (v-s)/s)
        #  as well as not take off the edges
        i = i + 1
    if verbose:
        for i in range(len(feeds)):
            print(f"fdnum,Tsys   {feeds[i]:2d}  {tsys[i]:10.5f}")
        print(f"<Tsys>  {np.nanmean(tsys):.5f} +/- {np.nanstd(tsys):.5f}")
        print(f"mode={mode}")
        print("TCAL=", tcal)

    return tsys


def plot_vegas(sdf, scans, title=None, tsys=False, inverse=False, edge=50, ylim=None):
    """
    plot vegas 16 beams

    Parameters
    ----------
    sdf : `GBTFITSLoad`
        data handle, containing one or more SDFITS files specific to GBT
    scans : list of ints
        DESCRIPTION.
    title : string, optional
        DESCRIPTION. The default is None.
    tsys : boolean, optional
        DESCRIPTION. The default is False.
    inverse : boolean, optional
        DESCRIPTION. The default is False.
    edge : int, optional
        DESCRIPTION. The default is 50.
    ylim : list of two floats, optional
        DESCRIPTION. The default is None, which autoscales.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(4, 4, sharex="col", sharey="row", gridspec_kw={"hspace": 0, "wspace": 0})

    for r in range(4):
        for c in range(4):
            p = r * 4 + c
            ax[r, c].plot()
            # ax[r,c].set_xlabel(f"{r} {c} {p}")
            for s in scans:
                v1 = sdf.gettp(scan=s, fdnum=p, calibrate=True, cal=False).timeaverage()
                if tsys:
                    s1 = sdf.gettp(scan=s + 1, fdnum=p, calibrate=True, cal=False).timeaverage()
                    vs = s1.data / (v1.data - s1.data)
                    if inverse:
                        vs = 1 / vs
                else:
                    vs = v1.data
                ax[r, c].plot(vs[edge:-edge], label=f"{p}")
                if tsys:
                    # note wwe haven't put a proper channel / freq axis
                    nc = len(vs[edge:-edge])
                    fix = np.ones(nc)
                    ax[r, c].plot(fix, color="black")
                # ax[r,c].scatter([0,900],[vs[edge],vs[-edge]], label=f"{p}")

                if ylim is not None:
                    ax[r, c].set_ylim(ylim)
            ax[r, c].legend()
    if title is None:
        plt.suptitle(sdf.filenames()[0])
    else:
        plt.suptitle(title)
    plt.tight_layout()

    if tsys:
        print(f"Showing sky/(vane-sky) for scans={scans},scans+1")
    else:
        print(f"Showing total power for scans={scans}")


def getnod(sdf, scans, beams, ifnum=0, plnum=0, tsys=None):
    """
    fake getnod() based on alternating gettp() with averaging done internally
    use the real sdf.getnod() for final analysis
    sdf:   the sdfits handle
    scans: list of two scans for the nodding
    beams: list of two beams for the nodding
    ifnum: the ifnum to use
    plnum: the plnum to use
    Returns the two nodding spectra, caller is responsible for averaging them, e.g.
         sp1.average(sp2)

    Parameters
    ----------
    sdf : GBTFITSLoad`
        data handle, containing one or more SDFITS files specific to GBT
    scans : list of 2 ints
        DESCRIPTION.
    beams : list of 2 ints
        DESCRIPTION.
    ifnum : int, optional
        DESCRIPTION. The default is 0.
    plnum : int, optional
        DESCRIPTION. The default is 0.
    tsys : float or list of two floats, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    sp1 : `Spectrum`
        DESCRIPTION.
    sp2 : `Spectrum`
        DESCRIPTION.

    """

    if tsys is None:
        tsys = np.array([1.0, 1.0])
    if np.isscalar(tsys):
        tsys = np.array([tsys, tsys])
    if len(tsys) == 1:
        tsys = np.array([tsys, tsys])  # because np.isscalar(np.array([1])) is False !

    ps1_on = sdf.gettp(scan=scans[0], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    ps1_off = sdf.gettp(
        scan=scans[1], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
    ).timeaverage()
    sp1 = (ps1_on - ps1_off) / ps1_off * tsys[0]

    ps2_on = sdf.gettp(scan=scans[1], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    ps2_off = sdf.gettp(
        scan=scans[0], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
    ).timeaverage()
    sp2 = (ps2_on - ps2_off) / ps2_off * tsys[1]

    sp1.meta["TSYS"] = tsys[0]
    sp2.meta["TSYS"] = tsys[1]

    return (sp1, sp2)
