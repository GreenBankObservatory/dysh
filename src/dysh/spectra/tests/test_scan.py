import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits

import dysh.util as util
from dysh.fits import gbtfitsload


class TestPSScan:
    def test_tsys(self, data_dir):
        """
        Test that `getps` results in the same system temperature as GBTIDL.
        """
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps = sdf.getps(scan=152, ifnum=0, plnum=0, fdnum=0)
        tsys = tps[0].tsys

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]
        hdu.close()

        assert (tsys - gbtidl_tsys)[0] == 0.0

    def test_compare_with_GBTIDL(self, data_dir):
        """ """
        # get filenames
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # psscan is a ScanList.
        psscan = sdf.getps(scan=152, plnum=0, ifnum=0, fdnum=0)
        assert len(psscan) == 1
        psscan.calibrate()
        # psscan_tavg is a spectrum.
        psscan_tavg = psscan.timeaverage(weights="tsys")

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        psscan_gbtidl = table["DATA"][0]
        hdu.close()

        # Compare.
        diff = psscan_tavg.flux.value - psscan_gbtidl
        assert abs(np.nanmedian(diff)) < 1e-9
        assert psscan_tavg.meta["EXPOSURE"] == table["EXPOSURE"][0]
        assert psscan_tavg.meta["TSYS"] == pytest.approx(table["TSYS"][0], 1e-14)

    def test_compare_with_GBTIDL_2(self, data_dir):
        """
        Test `getps` when working with multiple scans.
        """
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782"
        sdf_file = f"{data_path}/TGBT21A_501_11_NGC2782.raw.vegas.A.fits"
        gbtidl_file = f"{data_path}/TGBT21A_501_11_getps_scans_156-158_ifnum_0_plnum_0_timeaverage.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)
        ps_scans = sdf.getps(scan=[156, 158], ifnum=0, plnum=0, fdnum=0)
        ta = ps_scans.timeaverage()

        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        hdu.close()
        # changed from 1E-14 because I don't know how gbtidl calculated avg tsys
        assert ta.meta["TSYS"] == pytest.approx(table["TSYS"], rel=5e-6)
        assert ta.meta["EXPOSURE"] == table["EXPOSURE"]
        assert np.all(np.abs(table["DATA"][0] - ta.flux.value) < 3e-7)

    def test_ps_with_selection(self, data_dir):
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782"
        sdf_file = f"{data_path}/TGBT21A_501_11_NGC2782.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_scans1 = sdf.getps(scan=[156, 158], ifnum=0, plnum=0, fdnum=0)
        # because selection ANDs the selection rules, if the user selects [156,158]
        # and then we select [157,159] under the hood the AND will be exlusive and
        # the getps will fail. There is no easy way around this without a lot
        # of klugy code.  See issue #230
        sdf.select(scan=[156, 157, 158, 159])
        ps_scans2 = sdf.getps(ifnum=0, plnum=0, fdnum=0)
        assert len(ps_scans1) == 2
        assert len(ps_scans2) == 2
        assert np.all(ps_scans1[0]._calibrated == ps_scans2[0]._calibrated)
        assert np.all(ps_scans1[1]._calibrated == ps_scans2[1]._calibrated)

    @pytest.mark.skip(reason="We need to update this to work with multifits and ScanBlocks")
    def test_baseline_removal(self, data_dir):
        # get filenames
        data_path = f"{data_dir}/AGBT17A_404_01"
        sdf_file = f"{data_path}/AGBT17A_404_01_scan_19_prebaseline.fits"
        gbtidl_file = f"{data_path}/AGBT17A_404_01_scan_19_postbaseline.fits"
        gbtidl_modelfile = f"{data_path}/AGBT17A_404_01_scan_19_bline_model.fits"

        # generate the dysh result
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # has already been calibrated and smoothed slightly in gbtidl
        psscan = sdf.getspec(0)
        # 3rd order baseline fitted in gbtidl with regions [100,400] and [450,750]
        psscan.baseline(3, [(0, 100), (400, 450), (750, 819)], remove=True, model="chebyshev")

        # load gbtidl result
        # assuming the pre-baseline spectrum is still the same (can still test this if need be)
        hdu = fits.open(gbtidl_file)
        gbtidl_post = hdu[1].data["DATA"][0]
        hdu.close()
        hdu = fits.open(gbtidl_modelfile)
        gbtidl_bline_model = hdu[1].data["DATA"][0]  # noqa: F841
        hdu.close()

        diff = psscan - gbtidl_post  # noqa: F841

        # check that the spectra are the same but this won't pass right now
        # what is the tolerance for not passing?
        # [TODO] Find the right threshold for this
        # assert np.nanmean(diff) == 0

    def test_blank_integrations(self, data_dir):
        """
        Test `getps` when there are blanked integrations.
        """
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782_blanks"
        sdf_file = f"{data_path}/NGC2782.raw.vegas.A.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_sb = sdf.getps(scan=[156], fdnum=0, plnum=0, ifnum=0)
        ta1 = ps_sb.timeaverage()
        if False:
            # This should not raise any errors.
            ta2 = ps_sb.timeaverage(None)  # noqa: F841
            # Check if the time average is all NaNs.
            all_nan = np.isnan(ta1.flux.value).sum() == len(ta1.flux)
            assert ~all_nan
            # Check that the metadata is accurate.
            # The system temperature is different because of the squared averaging.
            assert abs(ps_sb[0].getspec(0).meta["TSYS"] - ta1.meta["TSYS"]) < 5e-16
            assert (ps_sb[0].getspec(0).meta["EXPOSURE"] - ta1.meta["EXPOSURE"]) == 0.0
            # Check if the time averaged data matches that from the first integration.
            # assert np.all(abs(ps_sb[0].getspec(0).flux.value - ta1[0].flux.value) < 2e-19)
            # Set to 5E-16 because Windows OS tests fail below that.  Need to understand why.
            assert np.all(abs(ps_sb[0].getspec(0).flux.value - ta1.flux.value) < 5e-16)

    def test_scan_write(self, data_dir, tmp_path):
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782_blanks"
        sdf_file = f"{data_path}/NGC2782.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_sb = sdf.getps(scan=[156], fdnum=0, plnum=0, ifnum=0)
        o = tmp_path / "scan_write"
        o.mkdir()
        testfile = o / "test_scan_write.fits"
        ps_sb[0].write(testfile, overwrite=True)

    def test_flag_write(self, data_dir, tmp_path):
        data_path = f"{data_dir}/AGBT18B_354_03"
        sdf_file = f"{data_path}/AGBT18B_354_03.raw.vegas"
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)
        sdf.clear_flags()
        flagarray = [60000, 90000]
        flagrange = np.arange(flagarray[0], flagarray[1] + 1)
        sdf.flag(scan=7, channel=[flagarray], intnum=[1])
        sdf.apply_flags()
        scan_block = sdf.getps(scan=6, ifnum=0, plnum=0, fdnum=0)
        output = tmp_path / "test_sb_flag_write.fits"
        # scanblock flags test
        scan_block.write(output, overwrite=True, flags=True)
        sdfin = gbtfitsload.GBTFITSLoad(output, flag_vegas=False)
        assert np.all(np.where(sdfin._sdf[0]._flagmask[0][1])[0] == flagrange)
        # scan flags test
        sb = sdfin.gettp(scan=7, ifnum=0, plnum=0, fdnum=0)
        assert np.all(np.where(sb[0]._calibrated.mask[1])[0] == flagrange)
        scanout = tmp_path / "test_scan_flag_write.fits"
        sb[0].write(scanout, overwrite=True, flags=True)
        scanin = gbtfitsload.GBTFITSLoad(scanout, flag_vegas=False)
        assert np.all(np.where(scanin._sdf[0]._flagmask[0][1])[0] == flagrange)
        tpscan = scanin.gettp(scan=7, ifnum=0, plnum=0, fdnum=0)
        assert np.all(np.where(tpscan[0]._calibrated.mask[1])[0] == flagrange)
        # try with flags=false
        scanout = tmp_path / "test_scan_flag_write2.fits"
        sb[0].write(scanout, overwrite=True, flags=False)
        scanin = gbtfitsload.GBTFITSLoad(scanout, flag_vegas=False)
        tpscan = scanin.gettp(scan=7, ifnum=0, plnum=0, fdnum=0)
        assert len(np.where(tpscan[0]._calibrated.mask[1])[0]) == 0

    def test_scale_units(self, data_dir):
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782_blanks"
        sdf_file = f"{data_path}/NGC2782.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tscale = "Flux"
        tunit = "Jy"
        ps_sb = sdf.getps(scan=156, fdnum=0, plnum=0, ifnum=0, units=tscale, zenith_opacity=0.08)
        assert ps_sb.tunit == tunit
        assert ps_sb.tscale == tscale
        assert ps_sb.tscale_fac == pytest.approx(0.55117614)
        assert ps_sb[0].tunit == tunit
        assert ps_sb[0].tscale == tscale

        ps_jy = ps_sb.timeaverage()
        assert ps_jy.meta["BUNIT"] == tunit
        assert ps_jy.meta["TUNIT7"] == tunit
        assert ps_jy.flux.unit.to_string() == tunit
        ps_jy_i = ps_sb[0].getspec(0)
        assert ps_jy_i.meta["BUNIT"] == tunit
        assert ps_jy_i.meta["TUNIT7"] == tunit
        assert ps_jy_i.flux.unit.to_string() == tunit

    def test_tcal(self, data_dir):
        """
        Test for getps with t_cal argument.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_org = sdf.getps(scan=152, plnum=0, ifnum=0, fdnum=0).timeaverage()
        ps_cal = sdf.getps(scan=152, plnum=0, ifnum=0, fdnum=0, t_cal=1.0).timeaverage()
        assert ps_cal.meta["TSYS"] == pytest.approx(ps_org.meta["TSYS"] / ps_org.meta["TCAL"])
        assert ps_cal.meta["TCAL"] == 1.0

    def test_ap_eff_surf_err(self, data_dir):
        """check that user-provided aperture efficiency and surface error are propogated correctly"""
        # assume that if it works for PSSCan it works for all Scan types???
        data_path = f"{data_dir}/AGBT18B_354_03"
        sdf_file = f"{data_path}/AGBT18B_354_03.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb1 = sdf.getps(scan=6, ifnum=0, plnum=0, fdnum=0, surface_error=0 * u.micron, zenith_opacity=0, units="ta*")
        # tscale_fac time loss efficiency = 0.99 should be 1/aperture efficiency
        assert np.all(sb1[0].ap_eff * sb1.tscale_fac * 0.99) == 1
        sb2 = sdf.getps(scan=6, ifnum=0, plnum=0, fdnum=0, ap_eff=1, zenith_opacity=0, units="ta*")
        # tscale fac times loss efficiency=0.99 should give 1 in this case
        assert np.all(sb2.tscale_fac * 0.99 == 1)
        x = sb1.timeaverage()
        assert pytest.approx(x.meta["AP_EFF"] / 0.709684575408695) == 1
        assert x.meta["SURF_ERR"] == 0.0
        assert x.meta["SE_UNIT"] == "micron"
        sb1 = sdf.getps(scan=6, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.1, units="ta*")
        x = sb1.timeaverage()
        assert pytest.approx(x.meta["TAU_Z"] / 0.1) == 1

    def test_vane(self, data_dir):
        """Test for getps with vane."""
        data_path = f"{data_dir}/TGBT24B_615_01/TGBT24B_615_01.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(data_path)
        pssb = sdf.getps(scan=86, ifnum=0, plnum=0, fdnum=10, vane=84, zenith_opacity=0.14, t_atm=268.85)
        psta = pssb.timeaverage()
        stats = psta.stats()
        assert stats["mean"].value == pytest.approx(-1.1362548)
        assert stats["median"].value == pytest.approx(-1.13077666)
        assert stats["rms"].value == pytest.approx(0.46247303)
        assert psta.meta["TSYS"] == pytest.approx(362.06692565749125)
        assert psta.meta["EXPOSURE"] == 29.732832173595035
        assert psta.meta["DURATION"] > psta.meta["EXPOSURE"]
        assert psta.meta["TSCALE"] == "Ta*"
        assert psta.meta["TSCALFAC"] == 1.0
        psta = sdf.getps(
            scan=86, ifnum=0, plnum=0, fdnum=10, vane=84, zenith_opacity=0.14, t_atm=268.85, units="flux"
        ).timeaverage()
        assert psta.meta["TSCALE"] == "flux"

    def test_getsigref_vane(self, data_dir):
        """Test for getsigref with vane."""
        data_path = f"{data_dir}/TGBT24B_615_01/TGBT24B_615_01.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(data_path)
        srsb = sdf.getsigref(scan=86, ifnum=0, plnum=0, fdnum=10, vane=84, zenith_opacity=0.14, t_atm=268.85, ref=87)
        srta = srsb.timeaverage()
        stats = srta.stats()
        result = {
            "mean": -1.14103931,
            "median": -1.13494246,
            "rms": 0.46420984,
            "min": -7.71431927,
            "max": 10.70582089,
        }
        for k in result:
            assert stats[k].value == pytest.approx(result[k])
        assert stats["nan"] == 31
        assert srta.meta["TSYS"] == pytest.approx(362.0670888626437)
        assert srta.meta["EXPOSURE"] == 58.5011952833595
        assert srta.meta["DURATION"] > srta.meta["EXPOSURE"]


class TestSubBeamNod:
    def test_compare_with_GBTIDL(self, data_dir):
        # get filenames
        # We still need a data file with a single scan in it
        sdf_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"
        gbtidl_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_snodka_43_ifnum_0_plnum_0_fdnum_1.fits"

        # Generate the dysh result.
        # snodka-style. Need test for method='cycle'
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sbn = sdf.subbeamnod(
            scan=43,
            ifnum=0,
            fdnum=1,
            plnum=1,
            calibrate=True,
            weights="tsys",
            method="scan",
        )

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        nodscan_gbtidl = hdu[1].data["DATA"][0]
        hdu.close()

        # Compare.
        ratio = sbn[0].getspec(0).flux.value / nodscan_gbtidl
        # kluge for now since there is a small wavy pattern in
        # the difference at the ~0.06 K level
        assert np.nanmedian(ratio) <= 0.998
        # make sure call to timeaverage functions
        sbn.timeaverage()

    def test_nodiode(self, data_dir):
        """
        Test for SubBeamNodScan without noise diodes.
        """
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        sdf_file = f"{data_dir}/AGBT17B_456_03/AGBT17B_456_03.raw.vegas.testtrim.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)

        # Cycle mode.
        sbn = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0).timeaverage()

        assert sbn.data.std() == pytest.approx(0.00222391)
        assert sbn.meta["EXPOSURE"] == 2.402706191448039
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]
        assert sbn.meta["SCAN"] == 20
        assert sbn.meta["TSYS"] == 1.0

        sbn = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, t_sys=105.0).timeaverage()

        assert sbn.meta["EXPOSURE"] == 2.402706191448039
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]
        assert sbn.meta["SCAN"] == 20
        assert sbn.meta["TSYS"] == pytest.approx(105.0)

        # Scan mode.
        sbn = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, method="scan").timeaverage()

        assert sbn.meta["SCAN"] == 20
        assert sbn.meta["TSYS"] == 1.0
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]

        sbn = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, method="scan", t_sys=100.0).timeaverage()

        assert sbn.meta["SCAN"] == 20
        assert sbn.meta["TSYS"] == pytest.approx(100.0)
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]

        # Equal weights.
        sbn_eq = sdf.subbeamnod(
            scan=20, ifnum=0, fdnum=10, plnum=0, method="scan", t_sys=100.0, weights=None
        ).timeaverage()

        assert (sbn.data - sbn_eq.data).sum() > 0.15
        assert sbn_eq.meta["SCAN"] == 20
        assert sbn_eq.meta["TSYS"] == pytest.approx(100.0)
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]

        # Smooth reference.
        sbn_smref = sdf.subbeamnod(
            scan=20, ifnum=0, fdnum=10, plnum=0, method="scan", t_sys=100.0, smoothref=10
        ).timeaverage()
        s = slice(2000, 6000)  # Clean channels.

        assert (sbn.data - sbn_smref.data)[s].sum() == pytest.approx(-5.375981637276874)
        assert sbn_smref.meta["SCAN"] == 20
        assert sbn_smref.meta["TSYS"] == pytest.approx(100.0)
        assert sbn_smref[s].stats()["rms"].value == pytest.approx(0.17582152178367458)
        assert sbn.meta["DURATION"] > sbn.meta["EXPOSURE"]

    def test_tcal(self, data_dir):
        """
        Test for SubBeamNodScan with t_cal argument.
        """
        sdf_file = f"{data_dir}/AGBT17B_456_03/AGBT17B_456_03.raw.vegas.testtrim.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)

        # Cycle mode.
        sbn_org = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0).timeaverage()
        sbn_cal = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, t_cal=1.0).timeaverage()
        assert sbn_cal.meta["TSYS"] == pytest.approx(sbn_org.meta["TSYS"] / sbn_org.meta["TCAL"])
        assert sbn_cal.meta["TCAL"] == 1.0

        # Scan mode.
        sbn_org = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, method="scan").timeaverage()
        sbn_cal = sdf.subbeamnod(scan=20, ifnum=0, fdnum=10, plnum=0, method="scan", t_cal=1.0).timeaverage()
        assert sbn_cal.meta["TSYS"] == pytest.approx(sbn_org.meta["TSYS"] / sbn_org.meta["TCAL"])
        assert sbn_cal.meta["TCAL"] == 1.0

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_synth_spectra(self, data_dir):
        """Test subbeamnod using synthithic spectra."""
        # Load a file with subbeamnod observations.
        sdf_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)

        # Generate fake data.
        def gauss(x, a, s, c):
            return a * np.exp(-((x - c) ** 2.0) / (2 * s**2.0))

        # Make a fake data set.
        gain = 1e8
        tsys = 100
        tcont = 10
        tcal = sdf["TCAL"][0]
        a = 10
        s = 5
        c = 800
        x = np.arange(0, sdf["DATA"].shape[1], 1)
        nchan = sdf["DATA"].shape[1]

        # Signal and reference integrations.
        sig_mask = sdf["SUBREF_STATE"].to_numpy() == -1
        ref_mask = sdf["SUBREF_STATE"].to_numpy() == 1
        n_sig = sig_mask.sum()
        n_ref = ref_mask.sum()

        # Generate astro signal.
        np.random.seed(0)  # Fix the noise.
        sig_mod = np.random.normal(loc=0, scale=0.1, size=(n_sig, nchan)) + gauss(x, a, s, c) + tcont
        ref_mod = np.random.normal(loc=0, scale=0.1, size=(n_ref, nchan))

        # Add noise diode.
        sig_cal_on = sdf["CAL"].to_numpy()[sig_mask] == "T"
        sig_mod[sig_cal_on] += tcal
        ref_cal_on = sdf["CAL"].to_numpy()[ref_mask] == "T"
        ref_mod[ref_cal_on] += tcal

        # Convert to counts.
        p_sig_mod = gain * (sig_mod + tsys)
        p_ref_mod = gain * (ref_mod + tsys)

        # Replace data.
        new_data = np.empty_like(sdf["DATA"])
        new_data[sig_mask] = p_sig_mod
        new_data[ref_mask] = p_ref_mod
        with pytest.warns(UserWarning):
            sdf["DATA"] = new_data
            sdf["TCAL"] = tcal

        # Calibrate.
        sbn_cycle = sdf.subbeamnod(scan=43, fdnum=1, plnum=1, ifnum=0).timeaverage()
        sbn_scan = sdf.subbeamnod(scan=43, fdnum=1, plnum=1, ifnum=0, method="scan").timeaverage()

        # Check data over a frequency interval.
        s_sbn = slice(30.0 * u.GHz, 30.5 * u.GHz)
        rms_cycle = sbn_cycle[s_sbn].flux.std()
        rms_scan = sbn_scan[s_sbn].flux.std()

        # Compare continuum level.
        assert pytest.approx(sbn_cycle.data.mean(), rms_cycle.value) == tcont
        assert pytest.approx(sbn_scan.data.mean(), rms_scan.value) == tcont

        # Compare exposure times.
        assert sbn_cycle.meta["EXPOSURE"] == 11.683568795384163
        assert sbn_scan.meta["EXPOSURE"] == 11.717374602730942

        # Compare system temperature.
        assert pytest.approx(sbn_scan.meta["TSYS"], rms_scan.value) == sbn_cycle.meta["TSYS"]
        assert pytest.approx(sbn_cycle.meta["TSYS"] - tcal / 2.0, rms_cycle.value) == tsys
        assert pytest.approx(sbn_scan.meta["TSYS"] - tcal / 2.0, rms_scan.value) == tsys

        # Compare RMS.
        assert pytest.approx(rms_cycle.value, abs=1e-2) == rms_scan.value

        # Number of rows.
        assert sdf.subbeamnod(scan=43, fdnum=1, plnum=1, ifnum=0, method="scan")[0].nrows == 96

        # Line amplitude.
        assert pytest.approx(sbn_cycle.data.max() - tcont, rms_cycle.value) == a
        assert pytest.approx(sbn_scan.data.max() - tcont, rms_scan.value) == a

    def test_vane(self, data_dir):
        """Test for subbeamnod with vane."""
        sdf_file = f"{data_dir}/AGBT18B_357_04/AGBT18B_357_04.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, skipflags=True)
        sbnsb = sdf.subbeamnod(scan=3, vane=1, ifnum=0, plnum=0, fdnum=10, zenith_opacity=0.1, t_atm=257.90)
        sbnta = sbnsb.timeaverage()
        stats = sbnta.stats()
        assert stats["mean"].value == pytest.approx(-0.09412889, abs=1e-4)
        assert stats["median"].value == pytest.approx(-0.08932108, abs=1e-4)
        assert stats["rms"].value == pytest.approx(0.43725554485804585, abs=1e-4)
        assert sbnta.meta["TSYS"] == pytest.approx(180.01411515504023, abs=1e-4)
        assert sbnta.meta["TSCALE"] == "Ta*"
        assert sbnta.meta["TSCALFAC"] == 1.0


class TestWeights:
    def test_weights(self, data_dir):
        sdf_file = f"{data_dir}/TGBT21A_501_11/NGC2782/TGBT21A_501_11_NGC2782.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb = sdf.getps(scan=[156, 158], ifnum=0, plnum=0, fdnum=0)
        x = sb.timeaverage(weights="tsys")
        assert np.all(x.weights - (sb._timeaveraged[0].weights + sb._timeaveraged[1].weights) == 0)
        sdf_file = f"{data_dir}/AGBT20B_014_03.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb = sdf.getfs(scan=6, fdnum=0, plnum=0, ifnum=0)
        w = np.ones((3, sb.nchan))
        x = sb.timeaverage(weights=w)
        assert x.weights[3071] == 2
        assert set(x.weights) == set([2.0, 3.0])

        x = sb.timeaverage(weights="tsys")
        assert np.ma.mean(x.weights) == pytest.approx(274.7940106210527)
        # now do a custom weight array
        scale = 100
        w = scale * np.random.rand(sb.nint, sb.nchan)
        x = sb.timeaverage(weights=w)
        # average weight should be roughly nint*np.mean(w) which should
        # also be 1/2 nint*scale.  1/2 because the mean of the random range (0,1) is 0.5
        assert np.mean(x.weights) == pytest.approx(sb.nint * np.mean(w), rel=1e-2)
        assert np.mean(x.weights) == pytest.approx(0.5 * sb.nint * scale, rel=1e-2)


class TestScanBase:
    def test_timeaverage_flags(self):
        """
        Test that `ScanBase.timeaverage()` produces the correct flags.
        """
        sdf_file = util.get_project_testdata() / "TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)
        channel = [2000, 6000]
        intnums = list(np.r_[0:70])
        chanslc = slice(channel[0], channel[1])
        sdf.flag(scan=152, channel=[channel], int=intnums)
        tp_sb = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0)
        # Mask is properly applied.
        assert np.all(tp_sb[0]._calibrated[intnums, chanslc].mask)
        tp = tp_sb[0].timeaverage()
        assert not np.all(tp[chanslc].mask)
        assert np.all(
            tp.weights[chanslc] == pytest.approx(tp_sb[0].tsys_weight.sum() - tp_sb[0].tsys_weight[intnums].sum())
        )
        assert np.all(
            tp[chanslc].weights == pytest.approx(tp_sb[0].tsys_weight.sum() - tp_sb[0].tsys_weight[intnums].sum())
        )
        assert np.all(tp.weights[: channel[0]] == pytest.approx(tp_sb[0].tsys_weight.sum()))
        # Channel selection in dysh is inclusive of the upper edge.
        assert np.all(tp.weights[channel[1] + 1 :] == pytest.approx(tp_sb[0].tsys_weight.sum()))


class TestTPScan:
    def test_len_and_units(self, data_dir):
        """
        Test that `TPScan` has the proper length and units.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)

        scan = 153
        tpsb = sdf.gettp(scan=153, ifnum=0, plnum=0, fdnum=0)

        assert len(tpsb[0]) == sdf.get_summary(scan=scan)["# INT"][0]
        assert tpsb.tscale == "Raw"
        assert tpsb.tunit == u.ct
        assert np.all(tpsb.tscale_fac == 1)

    def test_units_preserved(self, data_dir, tmp_path):
        """
        Test that TPScan preserves brightness units of calibrated data that was written out
        """
        data_path = f"{data_dir}/AGBT05B_047_01/AGBT05B_047_01.raw.acs"
        sdf_file = f"{data_path}/AGBT05B_047_01.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb = sdf.getps(scan=[51], ifnum=0, plnum=0, fdnum=0, units="flux", zenith_opacity=0.1)
        o = tmp_path / "tpsub"
        o.mkdir()
        testfile = o / "test_scanblock_write.fits"
        sb.write(fileobj=testfile, overwrite=True)
        sdf = gbtfitsload.GBTFITSLoad(testfile)
        tpsb = sdf.gettp(scan=[51], ifnum=0, plnum=0, fdnum=0)
        assert tpsb.tunit == sb.tunit
        assert tpsb.tscale == sb.tscale
        assert np.all(tpsb.tscale_fac == sb.tscale_fac)

    def test_tsys(self, data_dir):
        """
        Test that `gettp` produces the same system temperature as GBTDIL.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps = sdf.gettp(scan=153, ifnum=0, plnum=0, fdnum=0)
        tsys = tps[0].tsys

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]
        hdu.close()

        assert (tsys - gbtidl_tsys)[0] == 0.0

    def test_exposure(self, data_dir):
        """
        Test that `gettp` holds the same exposure times as GBTIDL.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tp = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0)

        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        hdu.close()

        for i in range(len(tp[0]._scanrows) // 2):
            assert tp[0].total_power(i).meta["EXPOSURE"] == table["EXPOSURE"][i]

    def test_compare_with_GBTIDL_tsys_weights(self, data_dir):
        """
        This test compares `gettp` when using radiometer weights.
        It takes a scan with multiple integrations and averages
        them using weights from the radiometer equation.
        It checks that the averaged exposure and spectra are the same,
        and that the system temperature is the same up to the precision
        used by GBTIDL.
        It also checks that the spectra, exposures and system temperatures are the
        same for individual integrations after calibrating them.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # tp is a ScanList
        tp = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0)
        assert len(tp) == 1
        # tpavg is a  Spectrum
        tpavg = tp.timeaverage()
        # assert len(tpavg) == 1

        # Check that we know how to add.
        assert tpavg.meta["EXPOSURE"] == tp[0].exposure.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]
        hdu.close()

        # Check exposure times.
        # The last row in the GBTIDL file is the averaged data.
        assert np.sum(table["EXPOSURE"][:-1] - tp[0].exposure) == 0.0
        assert tpavg.meta["EXPOSURE"] == table["EXPOSURE"][-1]
        # System temperature.
        # For some reason, the last integration comes out with a
        # difference in TSYS ~4e-10 rather than ~1e-14. Check why.
        assert np.all(abs(table["TSYS"][:-1] - tp[0].tsys) < 1e-9)
        assert abs(tpavg.meta["TSYS"] - table["TSYS"][-1]) < 1e-10
        # Data, which uses float -- 32 bits.
        assert np.sum(tp[0]._calibrated - data[:-1]) == 0.0
        assert np.nanmean((tpavg.flux.value - data[-1]) / data[-1].mean()) < 2**-32

    def test_compare_with_GBTIDL_equal_weights(self, data_dir):
        """
        This test compares `gettp` when using equal weights.
        It takes a scan with multiple integrations and averages
        them using unity weights.
        It checks that the exposure and averaged spectra are the same,
        and that the system temperature is the same up to the precision
        used by GBTIDL.
        """
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)
        tp = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0)
        assert len(tp) == 1
        # tpavg is a Spectrum
        tpavg = tp.timeaverage(weights=None)
        # assert len(tpavg) == 1

        # Check that we know how to add.
        assert tpavg.meta["EXPOSURE"] == tp[0].exposure.sum()
        assert tpavg.meta["DURATION"] == tp[0].duration.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]
        hdu.close()

        # Compare Dysh and GBTIDL.
        assert table["EXPOSURE"][0] == tpavg.meta["EXPOSURE"]
        assert table["DURATION"][0] == pytest.approx(tpavg.meta["DURATION"])
        assert abs(table["TSYS"][0] - tpavg.meta["TSYS"]) < 2**-32
        assert np.all((data[0] - tpavg.flux.value.astype(np.float32)) == 0.0)

    def test_t_sys_arg(self, data_dir):
        """
        Test for gettp with t_sys arg.
        """
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        t_sys = 10.0
        tpsb = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0, t_sys=t_sys)
        assert np.all(tpsb[0].tsys == 10.0)

    def test_t_cal_arg(self, data_dir):
        """
        Test for gettp with t_cal arg.
        """
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        t_cal = 10.0
        tpsb = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0, t_cal=t_cal)
        assert np.all(tpsb[0]._tcal == t_cal)
        tpsb_org = sdf.gettp(scan=152, ifnum=0, plnum=0, fdnum=0)
        assert np.all(tpsb_org[0].tsys / tpsb_org[0]._tcal * t_cal == tpsb[0].tsys)


class TestFSScan:
    def test_getfs_with_args(self, data_dir):
        sdf_file = f"{data_dir}/TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_504_01/TGBT21A_504_01.cal.vegas.fits"
        gbtidl_file_nofold = f"{data_dir}/TGBT21A_504_01/TGBT21A_504_01.nofold.vegas.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)

        fsscan = sdf.getfs(scan=20, ifnum=0, plnum=1, fdnum=0, fold=False)
        ta = fsscan.timeaverage(weights="tsys")
        #    we're using astropy access here, and use gbtfitsload.GBTFITSLoad() in the other test
        hdu = fits.open(gbtidl_file_nofold)
        table = hdu[1].data
        data = table["DATA"]
        hdu.close()
        level = 1e-6
        nm = np.nanmean(data[0] - ta.flux.value.astype(np.float32))
        assert abs(nm) <= level

        fsscan = sdf.getfs(scan=20, ifnum=0, plnum=1, fdnum=0, fold=True)
        ta = fsscan.timeaverage(weights="tsys")
        # We will be using GBTFITSLoad() here instead of astropy.
        if True:
            sdf2 = gbtfitsload.GBTFITSLoad(gbtidl_file)
            sp = sdf2.getspec(1).flux.value.astype(np.float32)
        else:
            hdu = fits.open(gbtidl_file)
            table = hdu[1].data
            data = table["DATA"]
            hdu.close()
            sp = data[1]
        level = 5e-3
        diff1 = sp - ta.flux.value.astype(np.float32)
        nm = np.nanmean(diff1[15000:20000])  # Use channel range around the line.
        assert abs(nm) <= level

        # Using interpolation to shift the data.
        fsscan = sdf.getfs(scan=20, ifnum=0, plnum=1, fdnum=0, fold=True, shift_method="interpolate")
        ta = fsscan.timeaverage(weights="tsys")
        diff2 = sp - ta.flux.value.astype(np.float32)
        nm = np.nanmean(diff2[15000:20000])
        assert abs(nm) <= level

        # Test with reference smoothing.
        fs_sb = sdf.getfs(scan=20, ifnum=0, plnum=0, fdnum=0, fold=True, smoothref=256)
        fs = fs_sb.timeaverage()
        assert fs.meta["EXPOSURE"] == pytest.approx(55.77325632268)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.meta["TSYS"] == pytest.approx(26.83285745447353)
        assert fs.stats()["mean"].value == pytest.approx(0.19299398411039015)
        assert fs.stats()["rms"].value == pytest.approx(5.4871938402480795)

        # nocal.
        fs_sb = sdf.getfs(scan=20, ifnum=0, plnum=0, fdnum=0, fold=True, smoothref=256, nocal=True)
        fs = fs_sb.timeaverage()
        assert fs.meta["EXPOSURE"] == pytest.approx(27.363056179345364)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.meta["TSYS"] == pytest.approx(1.0)
        assert fs.stats()["mean"].value == pytest.approx(0.007543638921500274)
        assert fs.stats()["rms"].value == pytest.approx(0.20901151397310902)

        # t_sys.
        t_sys = 124.0
        fs_sb = sdf.getfs(scan=20, ifnum=0, plnum=0, fdnum=0, fold=True, smoothref=256, t_sys=t_sys)
        fs = fs_sb.timeaverage()
        assert fs.meta["EXPOSURE"] == pytest.approx(55.77325632268)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.meta["TSYS"] == pytest.approx(t_sys)
        assert fs.stats()["mean"].value == pytest.approx(0.878913979866379)
        assert fs.stats()["rms"].value == pytest.approx(25.31681410111804)

        # nocal and t_sys.
        t_sys = 120.0
        fs_sb = sdf.getfs(scan=20, ifnum=0, plnum=0, fdnum=0, fold=True, smoothref=256, nocal=True, t_sys=t_sys)
        fs = fs_sb.timeaverage()
        assert fs.meta["EXPOSURE"] == pytest.approx(27.363056179345364)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.meta["TSYS"] == pytest.approx(t_sys)
        assert fs.stats()["mean"].value == pytest.approx(0.9052366705561161)
        assert fs.stats()["rms"].value == pytest.approx(25.081381667649197)

    def test_getfs_nocal(self):
        """
        Test for getfs without noise diode.
        """
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        sdf_file = util.get_project_testdata() / "AGBT20B_295_02/AGBT20B_295_02.raw.vegas.testtrim.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)

        # Test without system temperature.
        fs_sb = sdf.getfs(scan=12, ifnum=0, plnum=0, fdnum=10)
        assert fs_sb[0]._nocal
        fs = fs_sb.timeaverage()
        assert fs.meta["TSYS"] == 1.0
        assert fs.meta["EXPOSURE"] == pytest.approx(1.0926235028020896)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.stats()["mean"].value == pytest.approx(0.0011403941433353717)
        assert fs.stats()["rms"].value == pytest.approx(0.011686314570904896)

        # Test with system temperature.
        t_sys = 205.0
        fs_sb = sdf.getfs(scan=12, ifnum=0, plnum=0, fdnum=10, t_sys=t_sys)
        assert fs_sb[0]._nocal
        fs = fs_sb.timeaverage()
        assert fs.meta["TSYS"] == pytest.approx(t_sys)
        assert fs.meta["EXPOSURE"] == pytest.approx(1.0926235028020896)
        assert fs.stats()["mean"].value == pytest.approx(0.23378079938375124)
        assert fs.stats()["rms"].value == pytest.approx(2.395694487035504)

        # Test with reference smoothing.
        fs_sb = sdf.getfs(scan=12, ifnum=0, plnum=0, fdnum=10, smoothref=256)
        assert fs_sb[0]._nocal
        fs = fs_sb.timeaverage()
        assert fs.meta["TSYS"] == 1.0
        assert fs.meta["EXPOSURE"] == pytest.approx(2.115174908755242)
        assert fs.meta["DURATION"] > fs.meta["EXPOSURE"]
        assert fs.stats()["mean"].value == pytest.approx(0.0007369155360709178)
        assert fs.stats()["rms"].value == pytest.approx(0.010335320172352901)

    def test_tcal(self):
        """
        Test for getfs with t_cal argument.
        """
        sdf_file = (
            util.get_project_testdata() / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
        )
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        fs_org = sdf.getfs(scan=20, ifnum=0, plnum=1, fdnum=0).timeaverage()
        fs_cal = sdf.getfs(scan=20, ifnum=0, plnum=1, fdnum=0, t_cal=1.0).timeaverage()
        assert fs_cal.meta["TSYS"] == pytest.approx(fs_org.meta["TSYS"] / fs_org.meta["TCAL"])
        assert fs_cal.meta["TCAL"] == 1.0

    def test_vane(self):
        """
        Test for getfs with vane argument.
        """
        sdf_file = util.get_project_testdata() / "TGBT22A_603_05/TGBT22A_603_05.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        fs_sb = sdf.getfs(scan=12, ifnum=0, plnum=0, fdnum=2, vane=10, t_atm=273, zenith_opacity=0.1)
        assert fs_sb.tscale == "Ta*"
        fs = fs_sb.timeaverage()
        assert fs.meta["TSCALE"] == "Ta*"

    def test_mask_fold(self):
        """
        Test that FSScan handles masks.
        """
        sdf_file = (
            util.get_project_testdata() / "AGBT25A_504_03/AGBT25A_504_03.raw.vegas/AGBT25A_504_03.raw.vegas.A.fits"
        )
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)

        # Using interpolation for the shift.
        ta = sdf.getfs(scan=18, ifnum=0, plnum=1, fdnum=0, shift_method="interpolate").timeaverage()
        stats = ta.stats()
        assert stats["nan"] == 93
        assert np.isnan(ta.flux[5120].value)  # Single channel from the signal state.
        assert np.all(np.isnan(ta.flux[218:220].value))  # Two channels from the reference state.

        # Now using FFT shift.
        ta = sdf.getfs(scan=18, ifnum=0, plnum=1, fdnum=0, shift_method="fft").timeaverage()
        stats = ta.stats()
        assert stats["nan"] == 93
        assert np.isnan(ta.flux[5120].value)  # Single channel from the signal state.
        assert np.all(np.isnan(ta.flux[218:220].value))  # Two channels from the reference state.


class TestNodScan:
    def test_nodscan(self):
        """
        Test for `getnod` using data with noise diode.
        """
        # Reduce with dysh.
        # The SDFITS files used here did not flag vegas spurs, so don't flag them here
        fits_path = util.get_project_testdata() / "TGBT22A_503_02/TGBT22A_503_02.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path, flag_vegas=False)
        nod = sdf.getnod(scan=62, ifnum=0, plnum=0)
        assert len(nod) == 2
        nod_sp = nod.timeaverage()
        stats = nod_sp[int(2**15 * 0.1) : int(2**15 * 0.9)].stats()
        assert stats["rms"].value == pytest.approx(0.3502195954575188)
        assert stats["mean"].value == pytest.approx(0.21385395659416562)
        assert nod_sp.meta["DURATION"] > nod_sp.meta["EXPOSURE"]

    def test_tcal(self):
        """
        Test for getnod with t_cal argument.
        """
        fits_path = util.get_project_testdata() / "TGBT22A_503_02/TGBT22A_503_02.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        nod_sb_org = sdf.getnod(scan=62, ifnum=0, plnum=0)
        nod_sb_cal = sdf.getnod(scan=62, ifnum=0, plnum=0, t_cal=1.0)
        nod_org = nod_sb_org[0].timeaverage()
        nod_cal = nod_sb_cal[0].timeaverage()
        assert nod_cal.meta["TSYS"] == pytest.approx(nod_org.meta["TSYS"] / nod_org.meta["TCAL"])
        assert nod_cal.meta["TCAL"] == 1.0
        nod_org = nod_sb_org[1].timeaverage()
        nod_cal = nod_sb_cal[1].timeaverage()
        assert nod_cal.meta["TSYS"] == pytest.approx(nod_org.meta["TSYS"] / nod_org.meta["TCAL"])
        assert nod_cal.meta["TCAL"] == 1.0

    def test_vane(self):
        fits_path = util.get_project_testdata() / "AGBT22A_325_23/AGBT22A_325_23.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        nodsb = sdf.getnod(ifnum=0, plnum=0, vane=43, t_atm=265.48, zenith_opacity=0.21)
        assert nodsb[0].tsys.mean() == 195.85427050034397
        assert nodsb[1].tsys.mean() == 185.7013266638696
        assert nodsb.timeaverage().meta["TSCALE"] == "Ta*"
        assert nodsb.timeaverage().meta["TSCALFAC"] == 1.0
        nodsb = sdf.getnod(ifnum=0, plnum=0, vane=43, t_atm=265.48, zenith_opacity=0.21, units="flux")
        assert nodsb.timeaverage().meta["TSCALE"] == "flux"


class TestScanBlock:
    def test_scanblock_write_read(self, tmp_path):
        file = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/"
        g = gbtfitsload.GBTFITSLoad(file)
        sb = g.getps(scan=6, plnum=0, fdnum=0, ifnum=0)
        o = tmp_path / "sub"
        o.mkdir()
        testfile = o / "test_scanblock_write.fits"
        sb.write(fileobj=testfile, overwrite=True)
        g2 = gbtfitsload.GBTFITSLoad(testfile)
        x = g2.summary()  # simple check that basic function works.  # noqa: F841

    def test_baseline_subtraction(self, data_dir):
        """
        Tests for ScanBlock.subtract_baseline
        * Test that it works with a single Scan
        * Test that it works with multiple Scans
        * Test the tolerance limit
        * Test that undo_baseline works
        * Test that it raises an error with no calibrated data
        """
        data_path = f"{data_dir}/AGBT05B_047_01/AGBT05B_047_01.raw.acs"
        sdf_file = f"{data_path}/AGBT05B_047_01.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)

        # Single scan test.
        sb = sdf.getps(scan=[51], ifnum=0, plnum=0, fdnum=0)
        # check some properties while we're here
        assert sb.tunit == "K"
        assert sb.tscale == "Ta"
        assert np.all(sb.tscale_fac == 1)

        ta = sb.timeaverage()
        ta.baseline(exclude=[3000, 5000] * u.km / u.s, degree=1, remove=True)
        sb.subtract_baseline(ta.baseline_model)
        for s in sb:
            assert s.baseline_model == ta.baseline_model
            assert s.subtracted
        assert np.all(abs(ta.data - sb.timeaverage().data) < 1e-15)

        # Multiple scans.
        sb = sdf.getps(scan=[51, 53], ifnum=0, plnum=0, fdnum=0)
        ta = sb.timeaverage()
        ta.baseline(exclude=[3000, 5000] * u.km / u.s, degree=1, remove=True)

        # Check that tolerance limit works.
        with pytest.raises(ValueError):
            sb.subtract_baseline(ta.baseline_model, tol=0)
        sb.subtract_baseline(ta.baseline_model)
        for s in sb:
            assert s.baseline_model == ta.baseline_model
            assert s.subtracted
        # More integrations results in a larger difference wrt the
        # time average result, hence the change in threshold for the test.
        assert np.all(abs(ta.data - sb.timeaverage().data) < 1e-7)

        # Undo baseline for all scans.
        sb.undo_baseline()
        for s in sb:
            assert s.baseline_model is None
            assert not s.subtracted
        ta.undo_baseline()
        assert np.all(abs(ta.data - sb.timeaverage().data) < 1e-15)

        # No baseline allowed if not calibrated.
        sb = sdf.getps(scan=[51], ifnum=0, plnum=0, fdnum=0, calibrate=False)
        with pytest.raises(ValueError):
            sb.subtract_baseline(ta.baseline_model)

    def test_smooth(self, data_dir):
        data_path = f"{data_dir}/AGBT05B_047_01/AGBT05B_047_01.raw.acs"
        sdf_file = f"{data_path}/AGBT05B_047_01.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb = sdf.getps(scan=[51], ifnum=0, plnum=0, fdnum=0)
        rng = np.random.default_rng(12345)
        rdata = rng.random(sb[0]._calibrated.shape)
        rmask = sb[0]._calibrated.mask
        for width in [3, 5]:
            sb[0]._calibrated = np.ma.masked_array(rdata, rmask)

            mean = np.nanmean(sb[0]._calibrated)
            std = np.std(sb[0]._calibrated)
            sb[0].smooth(method="box", width=width, decimate=-1)
            hmean = np.nanmean(sb[0]._calibrated)
            hstd = np.std(sb[0]._calibrated)
            assert hmean == pytest.approx(mean, rel=1e-5)
            assert std / hstd == pytest.approx(np.sqrt(width), abs=1e-2)
            sb[0]._calibrated = np.ma.masked_array(rdata, rmask)
            sb[0].smooth(method="box", width=width, decimate=0)
            assert all(sb[0].delta_freq == np.array([x["CDELT1"] for x in sb[0].meta]))

    def test_tsys_weight(self, data_dir):
        data_path = f"{data_dir}/AGBT05B_047_01/AGBT05B_047_01.raw.acs"
        sdf_file = f"{data_path}/AGBT05B_047_01.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sb = sdf.getps(scan=[51], ifnum=0, plnum=0, fdnum=0)
        tsys_weights = sb.tsys_weight
        tsys_weights_scan = sb[0].tsys_weight
        assert all(tsys_weights[0] == tsys_weights_scan)

    def test_write_with_nan(self, data_dir, tmp_path):
        """Test that we can write a ScanBlock with NaN values."""
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file, flag_vegas=False)
        sb0 = sdf.getps(
            scan=152,
            fdnum=0,
            ifnum=0,
            plnum=0,
        )
        assert np.all(np.isnan(sb0[0]._get_all_meta("QD_EL")))
        with pytest.warns(UserWarning):
            sdf["QD_EL"] = 1
            sdf["QD_XEL"] = 1
        sb1 = sdf.getps(
            scan=152,
            fdnum=0,
            ifnum=0,
            plnum=0,
        )
        assert np.all(sb1[0]._get_all_meta("QD_EL") == [1])
        assert np.all(sb1[0]._get_all_meta("QD_XEL") == [1])
        sb0.append(sb1[0])
        sb0.write(tmp_path / "test.fits", overwrite=True)
