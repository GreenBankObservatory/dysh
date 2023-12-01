"""Load SDFITS files produced for GBO's 20m telescope """

from dysh.fits.sdfitsload import SDFITSLoad


class GB20MFITSLoad(SDFITSLoad):
    def __init__(self, filename, source=None, hdu=None, **kwargs):
        SDFITSLoad.__init__(self, filename, source, hdu)
        self._fix_columns()

    def _str_idx(self, text, bad_char="\x00"):
        """ """
        try:
            idx = text.index(bad_char)
        except ValueError:
            idx = None
        return idx

    def _fix_columns(self):
        """ """

        bad_char = "\x00"

        for i in range(len(self._bintable)):
            for j in range(self._binheader[i]["NAXIS2"]):
                text = self._bintable[i].data[j]["OBJECT"]
                idx = self._str_idx(text, bad_char="\x00")
                self._bintable[i].data[j]["OBJECT"] = text[:idx]
                text = self._bintable[i].data[j]["OBSERVER"]
                idx = self._str_idx(text, bad_char="\x00")
                self._bintable[i].data[j]["OBSERVER"] = text[:idx]
