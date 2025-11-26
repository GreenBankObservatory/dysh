import pathlib

import astropy.units as u
import pandas as pd
import pytest
from astropy.time import Time

import dysh
from dysh import util
from dysh.fits import gbtfitsload
from dysh.util.files import dysh_data
from dysh.util.selection import Selection, Flag

dysh_root = pathlib.Path(dysh.__file__).parent.resolve()


class TestSelectionWithPartialData:
    """Test Selection with partial data (e.g., from .index files)"""

    def test_selection_with_missing_alias_targets(self):
        """
        Test that Selection can be created when default alias targets are missing.

        Regression test for issue where loading from .index files failed because
        CRVAL1, CRVAL2, CRVAL3, SUBREF_STATE columns don't exist in the GBTIDL
        index file format.
        """
        # Simulate columns from a GBTIDL .index file (missing CRVAL1/2/3, SUBREF_STATE)
        df = pd.DataFrame(
            {
                "SCAN": [1, 2, 3],
                "OBJECT": ["test", "test2", "test3"],
                "ELEVATIO": [45.0, 50.0, 55.0],
                "PLNUM": [0, 0, 1],
                "DATE-OBS": ["2024-01-01T00:00:00", "2024-01-01T00:01:00", "2024-01-01T00:02:00"],
            }
        )

        # This should NOT raise an error
        sel = Selection(df)

        # Only aliases for existing columns should be set up
        assert "ELEVATION" in sel._aliases
        assert "SOURCE" in sel._aliases
        assert "POL" in sel._aliases
        assert sel._aliases["ELEVATION"] == "ELEVATIO"
        assert sel._aliases["SOURCE"] == "OBJECT"
        assert sel._aliases["POL"] == "PLNUM"

        # Aliases for missing columns should NOT be set up
        assert "FREQ" not in sel._aliases  # crval1 doesn't exist
        assert "RA" not in sel._aliases  # crval2 doesn't exist
        assert "DEC" not in sel._aliases  # crval3 doesn't exist
        assert "SUBREF" not in sel._aliases  # subref_state doesn't exist

    def test_flag_with_missing_alias_targets(self):
        """Test that Flag can also be created with missing alias targets."""
        df = pd.DataFrame(
            {
                "SCAN": [1, 2],
                "OBJECT": ["src1", "src2"],
                "ELEVATIO": [30.0, 35.0],
                "DATE-OBS": ["2024-01-01T00:00:00", "2024-01-01T00:01:00"],
            }
        )

        # This should NOT raise an error
        flag = Flag(df)

        # Only aliases for existing columns should be set up
        assert "ELEVATION" in flag._aliases
        assert "SOURCE" in flag._aliases

    def test_selection_works_with_existing_aliases(self):
        """Test that aliases that DO exist still work correctly."""
        df = pd.DataFrame(
            {
                "SCAN": [1, 2, 3, 4],
                "OBJECT": ["NGC1", "NGC1", "NGC2", "NGC2"],
                "ELEVATIO": [30.0, 35.0, 40.0, 45.0],
                "PLNUM": [0, 1, 0, 1],
                "DATE-OBS": ["2024-01-01T00:00:00", "2024-01-01T00:01:00", "2024-01-01T00:02:00", "2024-01-01T00:03:00"],
            }
        )

        sel = Selection(df)

        # Test using alias 'source' for 'OBJECT'
        sel.select(source="NGC1")
        assert len(sel.final) == 2

        sel.clear()

        # Test using alias 'pol' for 'PLNUM'
        sel.select(pol=0)
        assert len(sel.final) == 2


class TestSelection:
    """Test data selection"""

    def setup_method(self):
        self.file = util.get_project_testdata() / "TGBT21A_501_11/testselection.fits"
        self.sdf = gbtfitsload.GBTFITSLoad(self.file)
        self.gbtidl_flag_file = ""

    def test_selection_class(self):
        """
        Test that the Selection class selects as expected.
        We check that on each selection, the number of selected rows
        matches the known correct answer.
        """
        s = self.sdf.selection
        # polarization 0 and source name
        s.select(object="NGC2415", plnum=0)
        assert len(s._selection_rules[0]) == 5
        # Test that a duplicate selection results
        # in a warning.  Note here, we are
        # also testing that the alias 'pol' is interpreted
        # as plnum.
        with pytest.warns(UserWarning):
            s.select(object="NGC2415", pol=0, tag="this will warn", check=True)
        s.select(ifnum=[0, 2], tag="ifnums")
        assert len(s._selection_rules[1]) == 26
        # the AND of the selection rules becomes the final
        # selection
        assert len(s.final) == 3

        # test s.remove by both id and tag
        s.remove(0)
        assert len(s._selection_rules) == 1
        s.remove(tag="ifnums")
        assert len(s._selection_rules) == 0

        # This has the same final result than selections on separate lines
        s.select(object="NGC2415", plnum=0, ifnum=[0, 2])
        assert len(s.final) == 3
        s.clear()

        # Test using bool or char for CAL and SIG columns.
        s.select(cal="T")
        assert len(s.final) == 20
        assert set(s.final["CAL"]) == {"T"}
        s.clear()
        s.select(cal=True)
        assert len(s.final) == 20
        assert set(s.final["CAL"]) == {"T"}
        s.clear()
        s.select(sig="T")
        assert len(s.final) == 50
        s.clear()
        s.select(sig=True)
        assert len(s.final) == 50
        s.clear()
        with pytest.warns(UserWarning):  # There's no sig=F data.
            s.select(sig=False)
        assert len(s.final) == 0
        s.clear()
        with pytest.warns(UserWarning):  # There's no sig=F data.
            s.select(sig="F")
        assert len(s.final) == 0
        s.clear()

        # test select_range
        # lower limit.
        # RAs range from 53.24 to 247.04 degrees
        # Decs range from -15.64 to 54.57 degrees
        s.select_range(ra=(114,))  # default is degrees
        assert len(s.final) == 40
        # check that quantities work
        s.select_range(dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # again selection on the same line has a same final result
        s.clear()
        s.select_range(ra=(114,), dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # test select_within
        s.clear()
        # also verify that the selection variable name is
        # case insensitive
        # also note we aliased elevation for elevatio!
        s.select_within(eLEVaTIon=(18.0, 2))
        assert len(s.final) == 13

        # test selecting with a Time object.
        s.clear()
        t1 = Time("2021-02-10T08:00", scale="utc")
        t2 = Time("2021-02-10T09:00", scale="utc")
        s.select_range(utc=(t1, t2))
        assert len(s.final) == 10

        # test that an invalid  object for utc raise exception
        with pytest.raises(ValueError):
            s.select_range(utc=["asdad", 123])
        # np.datetime64 should also work
        s.clear()
        s.select_range(utc=(t1.datetime64, t2.datetime64))
        assert len(s.final) == 10
        # as should datetime
        s.clear()
        s.select_range(utc=(t1.datetime, t2.datetime))
        assert len(s.final) == 10
        # test select_channel
        a = [1, 4, (30, 40)]
        s.clear()
        s.select_channel(a)
        assert s._channel_selection == a
        assert len(s.final) == len(s)
        with pytest.raises(Exception):
            s.select_channel(["10", "a", 103])

    def test_flag_class(self):
        """
        Test that the Selection class selects as expected.
        We check that on each selection, the number of selected rows
        matches the known correct answer.
        """
        s = self.sdf.flags
        # polarization 0 and source name
        s.flag(object="NGC2415", plnum=0)
        assert len(s._selection_rules[0]) == 5
        # Test that a duplicate selection results
        # in a warning.  Note here, we are
        # also testing that the alias 'pol' is interpreted
        # as plnum.
        with pytest.warns(UserWarning):
            s.flag(object="NGC2415", pol=0, tag="this will warn", check=True)
        s.flag(ifnum=[0, 2], tag="ifnums")
        assert len(s._selection_rules[1]) == 26
        # the AND of the selection rules becomes the final
        # selection
        assert len(s.final) == 3

        # test s.remove by both id and tag
        s.remove(0)
        s.remove(tag="ifnums")
        assert len(s._selection_rules) == 0

        # This has the same final result than selections on separate lines
        s.flag(object="NGC2415", plnum=0, ifnum=[0, 2])
        assert len(s.final) == 3
        s.clear()

        # Test using bool or char for CAL and SIG columns.
        s.flag(cal="T")
        assert len(s.final) == 20
        assert set(s.final["CAL"]) == {"T"}
        s.clear()
        s.flag(cal=True)
        assert len(s.final) == 20
        assert set(s.final["CAL"]) == {"T"}
        s.clear()
        s.flag(sig="T")
        assert len(s.final) == 50
        s.clear()
        s.flag(sig=True)
        assert len(s.final) == 50
        s.clear()
        with pytest.warns(UserWarning):  # There's no sig=F data.
            s.flag(sig=False)
        assert len(s.final) == 0
        s.clear()
        with pytest.warns(UserWarning):  # There's no sig=F data.
            s.flag(sig="F")
        assert len(s.final) == 0
        s.clear()

        # test select_range
        # lower limit.
        # RAs range from 53.24 to 247.04 degrees
        # Decs range from -15.64 to 54.57 degrees
        s.flag_range(ra=(114,))  # default is degrees
        assert len(s.final) == 40
        # check that quantities work
        s.flag_range(dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # again selection on the same line has a same final result
        s.clear()
        s.flag_range(ra=(114,), dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # test select_within
        s.clear()
        # also verify that the selection variable name is
        # case insensitive
        # also note we aliased elevation for elevatio!
        s.flag_within(eLEVaTIon=(18.0, 2))
        assert len(s.final) == 13

        # test selecting with a Time object.
        s.clear()
        t1 = Time("2021-02-10T08:00", scale="utc")
        t2 = Time("2021-02-10T09:00", scale="utc")
        s.flag_range(utc=(t1, t2))
        assert len(s.final) == 10
        # test that an invalid  object for utc raise exception
        with pytest.raises(ValueError):
            s.flag_range(utc=["asdad", 123])
        # np.datetime64 should also work
        s.clear()
        s.flag_range(utc=(t1.datetime64, t2.datetime64))
        assert len(s.final) == 10
        # as should datetime
        s.clear()
        s.flag_range(utc=(t1.datetime, t2.datetime))
        assert len(s.final) == 10
        # test that a invalid object for utc raise exception
        with pytest.raises(ValueError):
            s.flag_range(utc=["asdad", 123])
        # test select_channel
        # a = [1, 4, (30, 40)]
        # s.clear()
        # s.select_channel(a)
        # assert s._channel_selection == a
        # assert len(s.final) == len(s)
        with pytest.raises(Exception):
            s.select_channel(["10", "a", 103])

    def test_flag_read(self):
        file = util.get_project_testdata() / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A6.fits"
        g = gbtfitsload.GBTFITSLoad(file, flag_vegas=False)
        assert len(g.flags._selection_rules) == 1
        assert g.flags._table["# SELECTED"][0] == 24
        assert g.flags._table["SCAN"][0] == "6"
        assert g.flags._table["IFNUM"][0] == "2"
        assert g.flags._table["FDNUM"][0] == "0"

    def test_flag_read_index_ok(self):
        # regression test for issue 457
        # first skipflags and read files separately
        f1 = dysh_data(test="AGBT22A_325_15/")
        sdf = gbtfitsload.GBTFITSLoad(f1, skipflags=True)
        flagA = f1 / "AGBT22A_325_15.raw.vegas.A.flag"
        flagB = f1 / "AGBT22A_325_15.raw.vegas.B.flag"
        sdf.flags.read(flagA)
        sdf.flags.read(flagB)
        # second, load flags at instantiation
        sdf = gbtfitsload.GBTFITSLoad(f1)
