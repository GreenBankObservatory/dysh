import os

from dysh.util import get_project_testdata
from dysh.util.download import from_url


class TestDownload:
    def setup_method(self):
        self.url = "https://github.com/GreenBankObservatory/dysh/raw/main/testdata/AGBT05B_047_01/gbtidl/AGBT05B_047_01.getps.acs.fits"

    def test_from_url(self):
        """Test that `from_url` is working"""

        try:
            fname = from_url(self.url)
            os.remove(fname)
        except Exception as exc:
            assert False, f"test_download_url raised the following exception: {exc}"

    def test_from_url_savepath(self):
        """Test that `from_url` returns the full path to the downloaded file."""

        path = get_project_testdata() / "AGBT05B_047_01/gbtidl"
        fname = from_url(self.url, path)

        assert fname == path / self.url.split("/")[-1]
