import numpy as np

from dysh.fits import gbtfitsload
from dysh.util.files import dysh_data


class TestVaneScan:
    def test_vane1(self):
        """
        test the vane/sky calibration of an Argus project
        using the VANECAL procedure
        """
        rtol = 1e-4  # do tests on this relative tolerance

        sdf_file = dysh_data(test="AGBT21B_024_14/AGBT21B_024_14_test")
        sdf2 = gbtfitsload.GBTFITSLoad(sdf_file)

        beam2 = sdf2.get_nod_beams(scan=331)
        assert len(beam2) == 2
        assert beam2[0] == 1 and beam2[1] == 9

        #  [1, 9]

        tcal = 272
        ifnum = 0
        plnum = 0
        tsys2 = sdf2.vanecal([329, 330], feeds=beam2, tcal=tcal, ifnum=ifnum, plnum=plnum)
        #  221.79946241 201.2739606

        sp1, sp2 = sdf2._getnod([331, 332], beam2, tsys=tsys2, ifnum=ifnum, plnum=plnum)
        sp3, sp4 = sdf2._getnod([333, 334], beam2, tsys=tsys2, ifnum=ifnum, plnum=plnum)
        sp5 = sp1.average([sp2, sp3, sp4])

        tint = sp5.meta["EXPOSURE"]
        assert np.isclose(tint, 3.942, rtol=rtol)  # 3.941774845123291

        tsys = sp5.meta["TSYS"]
        assert np.isclose(tsys, 210.54, rtol=rtol)  # 210.54325052318202

        rms = sp5.stats()["rms"].value  #  0.2357084722315336
        assert np.isclose(rms, 0.47509, rtol=rtol)

    def test_calseq2(self):
        """
        test the vane calibration of a W-band observation
        using the CALSEQ procedure
        """
        rtol = 1e-4  # do tests on this relative tolerance
        sdf_file = dysh_data(test="AGBT15B_244_07/AGBT15B_244_07_test")
        sdf3 = gbtfitsload.GBTFITSLoad(sdf_file)

        beam3 = sdf3.get_nod_beams(scan=131)
        assert len(beam3) == 2
        assert beam3[0] == 0 and beam3[1] == 1

        tsys3, g = sdf3.calseq(130, ifnum=1, plnum=0)  # 103.0501604820638
        assert np.isclose(tsys3, 104.33, rtol=rtol)

        sp1, sp2 = sdf3._getnod([131, 132], beam3, ifnum=1, plnum=0, tsys=tsys3)
        sp3 = sp1.average(sp2)

        tint = sp3.meta["EXPOSURE"]  # 11.939691305160522
        assert np.isclose(tint, 11.94, rtol=rtol)

        rms = sp3.stats()["rms"].value  # 0.2138536948704338
        assert np.isclose(rms, 0.21651, rtol=rtol)
