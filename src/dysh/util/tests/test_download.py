import os

from dysh.util.download import from_url


def test_download_url():
    """Test that wget is working"""

    url = "https://github.com/GreenBankObservatory/dysh/raw/main/testdata/AGBT05B_047_01/gbtidl/AGBT05B_047_01.getps.acs.fits"
    try:
        fname = from_url(url)
        os.remove(fname)
    except Exception as exc:
        assert False, f"test_download_url raised the following exception: {exc}"
