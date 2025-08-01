from pathlib import Path

import dysh.util.files as duf


class TestUtil:
    """Test dysh.util files functions"""

    def test_dysh_data(self):
        """Test dysh_data   (only limited testing possible)"""
        assert duf.dysh_data() == None  # noqa: E711
        # test=
        f1 = duf.dysh_data(test="getps")
        assert f1.exists() == True  # noqa: E712
        # example=
        #   given the dynamic nature, not sure if we should test example="getps"
        #   since it needs either $DYSH_DATA or /home/dysh/example_data
        # acccept=
        #   skipping
        # sdfits=
        #   skipping
        # dysh_data=
        assert duf.dysh_data("foo.fits", dysh_data="/tmp") == None  # noqa: E711
        #   this assume DYSH_DATA is not present
        f2 = duf.dysh_data(example="getps")
        assert f2.name == Path("AGBT05B_047_01.raw.acs.fits").name
