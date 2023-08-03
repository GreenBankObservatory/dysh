"""Load SDFITS files produced for GBO's 20m telescope """

from dysh.fits.sdfitsload import SDFITSLoad


class GB20MFITSLoad(SDFITSLoad):

    def __init__(self, filename, source=None, hdu=None, **kwargs):
        SDFITSLoad.__init__(self, filename, source, hdu)


    def _fix_object(self):
        """
        """

        for i in range(len(self._bintable)):
            for j in range(self._binheader[i]["NAXIS2"]):
                text = self._bintable[i].data[j]["OBJECT"]
                self._bintable[i].data[j]["OBJECT"] = text[:text.index("\x00")]
