import os

import wget


def test_flag_wget():
    """Test that wget is working"""

    url = "https://github.com/GreenBankObservatory/dysh/raw/main/testdata/AGBT05B_047_01/gbtidl/AGBT05B_047_01.getps.acs.fits"
    try:
        fname = wget.download(url)
        os.remove(fname)
    except Exception as exc:
        assert False, f"test_flag_wget raised the following exception: {exc}"
