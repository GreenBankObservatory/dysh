import numpy as np

import dysh.util.files as duf
from pathlib import Path

class TestUtil:
    """Test dysh.util files functions"""

    def test_dysh_data(self):
        """Test dysh_data   (only limited testing possible)"""
        assert duf.dysh_data() == None
        # test=
        f1 = duf.dysh_data(test="getps")
        assert f1.exists() == True
        # example=
        #   given the dynamic nature, not sure if we should test example="getps"
        #   since it needs either $DYSH_DATA or /home/dysh/example_data
        # acccept=
        #   skipping
        # sdfits=
        #   skipping
        # dysh_data=
        assert duf.dysh_data("foo.fits", dysh_data="/tmp") == None
        #   this assume DYSH_DATA is not present
        f2 = duf.dysh_data(example="test1")
        assert f2 == Path("AGBT05B_047_01.raw.acs.fits")
