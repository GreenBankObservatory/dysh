"""Load SDFITS files produced for GBO's 20m telescope"""

import numpy as np
from astropy import units as u

from dysh.coordinates import Observatory
from dysh.fits.sdfitsload import SDFITSLoad
from dysh.spectra import Spectrum


class GB20MFITSLoad(SDFITSLoad):
    def __init__(self, filename, source=None, hdu=None, **kwargs):
        SDFITSLoad.__init__(self, filename, source, hdu, index=False)
        self._fix_columns()
        self.create_index()
        self.selected_index = None
        self.selected_data = None

    def _str_idx(self, text, char="\x00"):
        """
        Find the index of `char` in `text`.
        """
        try:
            idx = text.index(char)
        except ValueError:
            idx = None
        return idx

    def _fix_columns(self):
        """
        Remove trailing characters from OBSMODE column and replaces
        OBJECT and OBSERVER calumns with the header values.
        """

        bad_char = b"\x00"

        for i in range(len(self._bintable)):
            for col in ["OBJECT", "OBSERVER"]:
                self._bintable[i].data[col] = self._hdu[i].header[col]
            for j in range(self._binheader[i]["NAXIS2"]):
                text = self._bintable[i].data[j]["OBSMODE"]
                idx = self._str_idx(text, char=bad_char)
                self._bintable[i].data[j]["OBSMODE"] = text[:idx]

    def _find_cal_idx(self):
        """ """

        if self.selected_index is None:
            raise TypeError("No data has been selected. Use GB20MFITSLoad.select()")

        swv_diff = np.diff(self.selected_index["SWPVALID"])
        jidxp = np.where(swv_diff == 1)[0]
        jidxm = np.where(swv_diff == -1)[0]
        self.cal_beg_idx = jidxp[0] + 1
        self.cal_end_idx = jidxm[1] + 1

    def _build_query(self, **kwargs):
        """ """

        query = ""
        for i, (k, v) in enumerate(kwargs.items()):
            kupp = k.upper()
            if kupp not in self._index.columns:
                continue
            query += f"{kupp} == {v}"
            if i < len(kwargs.items()) - 1:
                query += " and "

        return query

    def select(self, **kwargs):
        """ """

        query = self._build_query(**kwargs)

        self.selected_index = self._index.query(query)
        self.selected_data = self._bintable[0].data[list(self.selected_index.index)]
        self.selected_index.reset_index(inplace=True, drop=True)

    def _get_tsys_ps(self):
        """ """

        cal_beg_data = self.selected_data[: self.cal_beg_idx]
        tsys_beg = self._get_tsys(cal_beg_data)

        cal_end_data = self.selected_data[self.cal_end_idx :]
        tsys_end = self._get_tsys(cal_end_data)

        return np.nanmean([tsys_beg, tsys_end])

    def _get_tsys(self, cal_table):
        """ """

        on_mask = cal_table["CALSTATE"] == 1
        cal_on = np.nanmean(cal_table["DATA"][on_mask], axis=0)
        tcal = cal_table["TCAL"][on_mask].mean()
        of_mask = cal_table["CALSTATE"] == 0
        cal_of = np.nanmean(cal_table["DATA"][of_mask], axis=0)

        tsys = cal_of / (cal_on - cal_of) * tcal

        return np.nanmedian(tsys)

    def do_sigref(self, sig, ref, tsys):
        """ """

        return (sig - ref) / ref * tsys

    def _do_sigref_ps(self, tsys):
        """ """

        on_data = self.selected_data[list(self.selected_index.query("SWPVALID == 1 and OBSMODE == 'onoff:on'").index)]
        of_data = self.selected_data[list(self.selected_index.query("SWPVALID == 1 and OBSMODE == 'onoff:off'").index)]

        on_avg = np.average(on_data["DATA"], axis=0, weights=on_data["EXPOSURE"])
        of_avg = np.average(of_data["DATA"], axis=0, weights=of_data["EXPOSURE"])

        cal = self.do_sigref(on_avg, of_avg, tsys)

        return cal

    def make_meta_ps(self, tsys):
        """ """

        cols = self.selected_index.query("SWPVALID == 1 and OBSMODE == 'onoff:on'")
        meta = cols.iloc[0].to_dict()
        meta["TSYS"] = tsys
        meta["EXPOSURE"] = cols["EXPOSURE"].sum()
        meta["DURATION"] = cols["DURATION"].sum()
        meta["NAXIS1"] = int(meta["TDIM7"][1:-1].split(",")[0])
        meta["CUNIT1"] = "Hz"
        meta["CUNIT2"] = "deg"
        meta["CUNIT3"] = "deg"
        meta["CRVAL2"] = cols["CRVAL2"].mean()
        meta["CRVAL3"] = cols["CRVAL3"].mean()
        meta["RESTFRQ"] = meta["RESTFREQ"]

        return meta

    def getps(self, ifnum=0, plnum=0):
        r"""
        Calibrate position switched observations.
        It time averages the signal and reference spectra
        and then calibrates using:

        :math:`T_{\rm{sys}}\frac{\mathrm{SIG}}{\mathrm{SIG}-\mathrm{REF}}`

        Parameters
        ----------
        ifnum : int
            Spectral window to calibrate.
        plnum : int
            Polarization to calibrate.

        Returns
        -------
        cal : :class:`~spectra.spectrum.Spectrum`
            Calibrated spectrum.
        """

        self.select(ifnum=ifnum, plnum=plnum)
        self._find_cal_idx()
        tsys = self._get_tsys_ps()
        cal = self._do_sigref_ps(tsys)

        meta = self.make_meta_ps(tsys)

        return Spectrum.make_spectrum(cal * u.K, meta=meta, observer_location=Observatory["GB20M"])
