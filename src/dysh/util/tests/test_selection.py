import pathlib

import astropy.units as u
import pytest
from astropy.time import Time

import dysh
from dysh import util
from dysh.fits import gbtfitsload
from dysh.util.selection import Selection

dysh_root = pathlib.Path(dysh.__file__).parent.resolve()


class TestSelection:
    """Test data selection"""

    def setup_method(self):
        self.file = util.get_project_testdata() / "TGBT21A_501_11/testselection.fits"

    def test_selection_class(self):
        """
        Test that the Selection class selects as expected.
        We check that on each selection, the number of selected rows
        matches the known correct answer.
        """
        sdf = gbtfitsload.GBTFITSLoad(self.file)
        s = Selection(sdf)
        # the unfilter sdfits index should be the same as
        # the default selection except for addition of UTC column
        assert sdf._index.equals(s.drop(columns="UTC", inplace=False))
        # polarization 0 and source name
        s.select(object="NGC2415", plnum=0)
        assert len(s._selection_rules[0]) == 5
        # Test that a duplicate selection results
        # in a warning.  Note here, we are
        # also testing that the alias 'pol' is interpreted
        # as plnum.
        with pytest.warns(UserWarning):
            s.select(object="NGC2415", pol=0, tag="this will warn")
        s.select(ifnum=[0, 2], tag="ifnums")
        assert len(s._selection_rules[1]) == 26
        # the AND of the selection rules becomes the final
        # selection
        assert len(s.final) == 3

        # test s.remove by both id and tag
        s.remove(0)
        s.remove(tag="ifnums")
        assert len(s._selection_rules) == 0

        s = Selection(sdf)
        # This has the same final result than selections on separate lines
        s.select(object="NGC2415", plnum=0, ifnum=[0, 2])
        assert len(s.final) == 3

        # test select_range
        # lower limit.
        # RAs range from 53.24 to 247.04 degrees
        # Decs range from -15.64 to 54.57 degrees
        s = Selection(sdf)
        s.select_range(ra=(114,))  # default is degrees
        assert len(s.final) == 40
        # check that quantities work
        s.select_range(dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # again selection on the same line has a same final result
        s = Selection(sdf)
        s.select_range(ra=(114,), dec=[2400, 7500] * u.arcmin)
        assert len(s.final) == 20
        # test select_within
        s = Selection(sdf)
        # also verify that the selection variable name is
        # case insensitive
        # also note we aliased elevation for elevatio!
        s.select_within(eLEVaTIon=(18.0, 2))
        assert len(s.final) == 13

        # test selecting with a Time object.
        s = Selection(sdf)
        s.select_range(utc=(Time("2021-02-10T08:00", scale="utc"), Time("2021-02-10T09:00", scale="utc")))
        # test that a non-Time object for utc raise exception
        with pytest.raises(ValueError):
            s.select_range(utc=["asdad", 123])
        # test select_channel
        a = [1, 4, (30, 40)]
        s = Selection(sdf)
        s.select_channel(a)
        assert s._channel_selection == a
        assert len(s.final) == len(s)
        with pytest.raises(Exception):
            s.select_channel(["10", "a", 103])
