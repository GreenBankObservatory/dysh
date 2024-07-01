#!/usr/bin/env python
import copy
import sys

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter
from astropy.modeling.models import Gaussian1D
from astropy.modeling.polynomial import Polynomial1D
from astropy.table import Table
from astropy.units import cds
from astropy.wcs import WCS
from scipy.optimize import leastsq
from specutils import SpectralRegion, Spectrum1D, SpectrumList
from specutils.fitting import fit_continuum


def baseline(speclist, order, exclude=None, plot=False, maxspec=1000, scipy=False):
    # last = min(len(speclist),maxspec)
    print(f"BL {order} for {len(speclist)} spectra")
    # for p in speclist[0:last]:
    i = 0
    bad = 0
    model = Polynomial1D(degree=order)
    fitter = LinearLSQFitter(calc_uncertainties=True)
    # fitter = LevMarLSQFitter(calc_uncertainties=True)
    print(f"MODEL {model} FITTER {fitter}")
    try:
        if exclude is not None:
            for p in speclist:
                if np.isnan(p.data).all():
                    bad += 1
                    continue
                fc = fit_continuum(spectrum=p, model=model, fitter=fitter, exclude_regions=[exclude])
                i = i + 1
        else:
            for p in speclist:
                if np.isnan(p.data).all():
                    bad += 1
                    continue
                fc = fit_continuum(spectrum=p, model=model, fitter=fitter)
                i = i + 1
    except Exception as e:
        print(f"At spectrum {i}, Exception was {e}")
        print(p)
        print("DATA MEAN: ", np.nanmean(p.data))
        return p
    if plot:
        fig, ax = plt.subplots()
        ax.plot(x, p.flux)
        ax.plot(x, fc(x))
        plt.show()
    print(f"NUMBER OF BAD SPECTRA: {bad}")
    return None


def get_size(obj, seen=None):
    # https://goshippo.com/blog/measure-real-size-any-python-object/
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, "__dict__"):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def dcmeantsys(calon, caloff, tcal, mode=0, fedge=10, nedge=None):
    """
    following the GBTIDL routine with same name, get the tsys from
    the neighboring calon and caloff we define an extra way to set
    the edge size, nedge, if you prefer to use number of edge channels
    instead of the inverse fraction

    calon/caloff is meant to reflect the state of the noise diode

    mode=0     do the mean before the division
    mode=1     do the mean after the division
    """
    nchan = len(calon)
    if nedge is None:
        nedge = nchan // fedge  # 10 %
    if mode == 0:
        meanoff = np.mean(caloff[nedge:-nedge])
        meandiff = np.mean(calon[nedge:-nedge] - caloff[nedge:-nedge])
        meanTsys = meanoff / meandiff * tcal + tcal / 2.0
    else:
        meanTsys = np.mean(caloff[nedge:-nedge] / (calon[nedge:-nedge] - caloff[nedge:-nedge]))
        meanTsys = meanTsys * tcal + tcal / 2.0
    return meanTsys


# In[330]:


def uniq(seq):
    """from http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]


def sonoff(scan, procseqn):
    """
    return the list of On and Off scan numbers
    there must be a more elegant python way to do this....
    """
    sp = {}
    for i, j in zip(scan, procseqn):
        sp[i] = j

    us1 = uniq(scan)
    up1 = uniq(procseqn)

    sd = {}
    for i in up1:
        sd[i] = []

    for s in us1:
        sd[sp[s]].append(s)

    return sd


class Obsblock:
    """Class that holds a series of spectra on which bulk operations can be performed"""

    def __init__(self, speclist, index):
        self._speclist = speclist
        self._index = index  # pandas dataframe

    def __getitem__(self, i):
        return self._speclist[i]

    def __len__(self):
        return len(self._speclist)

    def __op__(self, opname):
        pass


# # SDFITSLoad
#
# This is the class that loads an SDFITS file. Normally not called by users, but by classes such as GBTLoadPS()
#
#

# In[334]:


class SDFITSLoad(object):
    """
    Container for a bintable(s) from selected HDU(s)
    """

    def __init__(self, filename, source=None, hdu=None, **kwargs):
        """ """
        print("==SDFITSLoad %s" % filename)
        cds.enable()  # to get mmHg
        kwargs_opts = {"fix": False}
        kwargs_opts = {"wcs": False}
        kwargs_opts.update(kwargs)
        self._filename = filename
        self._bintable = []
        self._ptable = []
        self._binheader = []
        self._data = []
        self._obsblock = []  # list of SpecList
        self._hdu = fits.open(filename)
        self._header = self._hdu[0].header
        self.load(hdu, **kwargs_opts)
        # self.load_pandas()
        # self._hdu.close()  # can't access hdu[i].data member of you do this.

    @property
    def filename(self):
        return self._filename

    def reset(self, hdu=None):
        self._bintable = []
        self._binheader = []
        self._ptable = []
        self._data = []
        self._spectra = []
        self._hdu = fits.open(self._filename)
        self._header = self._hdu[0].header
        # if hdu is not None:
        #   self.load(src,hdu)

    def load_pandas(self, hdu=None):
        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1, len(self._hdu))
        self._ptable = []
        for i in ldu:
            t = Table.read(self._hdu[i])
            t.remove_column("DATA")
            self.stripT(t)
            print(f"doing pandas for HDU {i}")
            self._ptable.append(t.to_pandas())
            del t

    def load(self, hdu=None, **kwargs):
        """
        for given hdu make this bintable available
        Note mmHg and UTC are unrecognized units.  mmHg is in astropy.units.cds but UTC is just wrong.
        """

        self._bintable = []
        self._ptable = []
        self._binheader = []
        self._data = []
        self._nrows = []
        source = kwargs.get("source", None)
        fix = kwargs.get("fix", False)
        wcs = kwargs.get("wcs", False)

        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1, len(self._hdu))
        for i in ldu:
            j = i - 1
            self._bintable.append(self._hdu[i])
            self._binheader.append(self._hdu[i].header)
            # TODO: don't allow preselection here, it just screws bookkeepingu p.
            # All selection must happen in obsblocks.
            # if source is None:
            # self._data.append(self._bintable[j].data[:]["DATA"])
            # else:
            #    wh1 = np.char.strip(self._bintable[j]['OBJECT']) == source # true/false array
            #    if wh1.sum() == 0: #all False none found
            #        srcs = np.unique(self._bintable[j]['OBJECT'])
            #        raise Exception(f"Source name {source} not found in HDU {i}. Sources present are {srcs}")
            #    self._data.append(self._bintable[j]["DATA"][wh1])
            #    mask = np.where(wh1)[0]
            #    self._bintable[j] = self._bintable[j][mask]  # header will be wrong?
            self._nrows.append(self._binheader[j]["NAXIS2"])

        if kwargs.get("index", False):
            self.load_pandas(hdu)

    def stripT(self, b):
        # remove leading and trailing chars from all strings in table
        for n in b.colnames:
            if np.issubdtype(b.dtype[n], str):
                b[n] = np.char.strip(b[n])

    def strip(self):
        # remove leading and trailing chars from all strings in table
        for b in self._ptable:
            for n in b.colnames:
                if np.issubdtype(b.dtype[n], str):
                    b[n] = np.char.strip(b[n])

    def _loadlists(self, hdu, fix=False, wcs=False, maxspect=1e16):
        self._obsblock = []
        i = 0
        k = -1
        print("HDU = ", hdu)
        if hdu is not None:
            b = self._ptable[hdu - 1]
            rawspect = self._bintable[i].data["DATA"]
            sl = SpectrumList()
            maxload = int(np.min([maxspect, self.nrows(i)]))
            print(f"Creating {maxload} Spectrum1D in bintable {i} HDU {hdu}", file=sys.stderr)
            for j in range(maxload):
                k = k + 1
                # need extra [[]] because we have 1x1 spatial NAXIS
                # otherwise, slicing the spectrum won't work.
                if wcs:
                    sp = np.array([[self.rawspectrum(i, j)]])
                else:
                    # sp = self.rawspectrum(i,j)*u.K
                    sp = np.copy(rawspect[j])  # *u.K
                    # sp = np.random.rand(32768)*u.K
                naxis1 = sp.shape[0]  # self.nchan(i)
                printme = int(0.1 * len(b))
                if (k % printme) == 0:
                    print(f"Row {k} nchan {naxis1} {type(sp)}", file=sys.stderr)
                    # print(f"NAXIS1 is {naxis1}",file=sys.stderr)
                crval1 = b["CRVAL1"][j]
                cdelt1 = b["CDELT1"][j]
                crpix1 = b["CRPIX1"][j]
                ctype1 = b["CTYPE1"][j]
                # Ensure rest frequency is in Hertz
                # CUNIT1 is not always present
                restfrq = b["RESTFREQ"][j]
                if "CUNIT1" in b.columns:
                    cunit1 = b["CUNIT1"][j]
                    rfq = restfrq * u.Unit(cunit1)
                    restfreq = rfq.to("Hz").value
                cunit1 = "Hz"
                crval2 = b["CRVAL2"][j]
                crval3 = b["CRVAL3"][j]
                ctype2 = b["CTYPE2"][j]
                ctype3 = b["CTYPE3"][j]
                # 'FREQ-OBS' to 'FREQ'; assuming SPECSYS='TOPOCENT'
                # if ctype1 == 'FREQ-OBS': ctype1  = 'FREQ'
                # only axis1 needs a full description, axis2,3,4 are all single points
                if wcs:
                    wcs = WCS(
                        header={
                            "CDELT1": cdelt1,
                            "CRVAL1": crval1,
                            "CUNIT1": cunit1,
                            "CTYPE1": "FREQ",
                            "CRPIX1": crpix1,
                            "RESTFRQ": restfrq,
                            "CTYPE2": ctype2,
                            "CRVAL2": crval2,
                            "CRPIX2": 1,
                            "CTYPE3": ctype3,
                            "CRVAL3": crval3,
                            "CRPIX3": 1,
                            "CUNIT2": "deg",
                            "CUNIT3": "deg",
                            "NAXIS1": naxis1,
                            "NAXIS2": 1,
                            "NAXIS3": 1,
                        },
                        fix=fix,
                    )
                else:
                    wcs = None
                # GBT really fucks up FREQ/VELDEF/VELFRAME
                # if False:
                if "VELFRAME" in b.columns:
                    vframe = b["VELFRAME"][j]
                elif "VFRAME" in b.columns:
                    vframe = b["VFRAME"][j]
                else:
                    vframe = None
                if "VELDEF" in b.columns:
                    vdef = b["VELDEF"][j]
                else:
                    vdef = None

                meta = {
                    "CDELT1": cdelt1,
                    "CRVAL1": crval1,
                    "CUNIT1": cunit1,
                    "CTYPE1": "FREQ",
                    "CRPIX1": crpix1,
                    "RESTFRQ": restfrq,
                    "CTYPE2": ctype2,
                    "CRVAL2": crval2,
                    "CRPIX2": 1,
                    "CTYPE3": ctype3,
                    "CRVAL3": crval3,
                    "CRPIX3": 1,
                    "CUNIT2": "deg",
                    "CUNIT3": "deg",
                    "NAXIS1": naxis1,
                    "NAXIS2": 1,
                    "NAXIS3": 1,
                    "VELDEF": vdef,
                    "VELFRAME": vframe,
                }

                convention = self.velocity_convention(vdef, vframe)
                # meta = dict(b.loc[j])# Necessary? Since we are sending whole pandas table to Obsblock
                if fix:
                    self.fix_meta(meta)
                if False:
                    try:
                        convention = self.velocity_convention(meta["VELDEF"], meta["VELFRAME"])
                    except Exception:
                        # print("WARNING: insufficient veldef/velframe, assuming convention is 'doppler_radio'")
                        convention = "doppler_radio"
                meta = {}
                sl.append(Spectrum1D(flux=sp * u.K, wcs=wcs, meta=meta, velocity_convention=convention))
            self._obsblock.append(Obsblock(sl, self._ptable[i]))
            i = i + 1

    def fix_meta(self, meta):
        """Do any repair to the meta/header for peculariaties in definitions from a particular observatory
        The passed-in dictionary will be repaired in place.
        At minimum this method must populate meta['VELDEF'] and meta['VELFRAME']
        """
        pass

    def velocity_convention(self, veldef, velframe):
        # sub-classes must implement this so I can keep this class generic.
        # GBT uses VELDEF and VELFRAME incorrectly.
        return "doppler_radio"

    def udata(self, bintable, key):
        return uniq(self._data[bintable][key])

    def ushow(self, bintable, key):
        print(f"{bintable} {key}: {self.udata(bintable,key)}")

    def naxis(self, bintable, naxis):
        nax = f"NAXIS{naxis}"
        return self._binheader[bintable][nax]

    def nintegrations(self, bintable, source=None):
        if source is not None:
            nint = np.shape(np.char.strip(self._data[bintable]["OBJECT"]) == source)[0] // self.npol(bintable)
        else:
            nint = np.shape(self._data[bintable])[0] // self.npol(bintable)
        return nint

    def rawspectra(self, bintable):
        return self._bintable[bintable].data[:]["DATA"]

    def rawspectrum(self, bintable, i):
        return self._bintable[bintable].data[:]["DATA"][i]

    def nrows(self, bintable):
        return self._nrows[bintable]

    def nchan(self, bintable):
        return np.shape(self.rawspectrum(bintable, 1))[0]

    def npol(self, bintable):
        return len(self.udata(bintable, "CRVAL4"))

    def sources(self, bintable):
        return self.udata(bintable, "OBJECT")

    def scans(self, bintable):
        return self.ushow(bintable, "SCAN")

    def __len__(self):
        return self.nrows

    def _summary(self, bintable):
        j = bintable
        nrows = self.naxis(j, 2)
        nflds = self._binheader[j]["TFIELDS"]
        restfreq = np.unique(self._data[j]["RESTFREQ"]) / 1.0e9
        #
        print("HDU       %d" % (j + 1))
        print("BINTABLE: %d rows x %d cols with %d chans" % (self._nrows[j], nflds, self.nchan(j)))
        print("Selected  %d/%d rows" % (self._nrows[j], nrows))
        print("Sources: ", self.sources(j))
        print("RESTFREQ:", restfreq, "GHz")
        print("Scans:   ", self.scans(j))
        print("Npol:    ", self.npol(j))
        print("Nint:    ", self.nintegrations(j))

    def summary(self):
        print("File:     %s" % self._filename)
        for i in range(len(self._bintable)):
            self._summary(i)

    def __repr__(self):
        return self._filename


class GBTLoad(SDFITSLoad):
    def __init__(self, filename, src=None, hdu=None):
        """
        Holds a raw "unstructured" series of scans, normally not used by users
        """
        SDFITSLoad.__init__(self, filename, src, hdu, fix=False)
        print("==GBTLoad %s" % filename)

        self.ushow(0, "OBJECT")
        self.ushow(0, "SCAN")
        self.ushow(0, "SAMPLER")
        # ushow('PLNUM')
        # ushow('IFNUM')
        self.ushow(0, "SIG")
        self.ushow(0, "CAL")
        self.ushow(0, "PROCSEQN")
        self.ushow(0, "PROCSIZE")
        self.ushow(0, "OBSMODE")
        self.ushow(0, "SIDEBAND")


if __name__ == "__main__":
    import sys
    GBTLoad(sys.argv[1])
