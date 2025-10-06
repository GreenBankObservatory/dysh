import matplotlib.pyplot as plt
import numpy as np


def plot_vegas(sdf, scans, title=None, tsys=False, inverse=False, edge=50, ylim=None):
    """
    plot the vegas 16 beams like in vanecal.pro
    This routine is still rudimentary for a quick-look.

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

    _fig, ax = plt.subplots(4, 4, sharex="col", sharey="row", gridspec_kw={"hspace": 0, "wspace": 0})

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
