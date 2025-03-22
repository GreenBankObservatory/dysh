"""Tests for dysh.config.core"""

from astropy.config import set_temp_config

import dysh.config as dc


class TestConfig:

    def test_get_config_dir_path(self):
        path = dc.get_config_dir_path(rootname="dysh-test")
        assert path.exists()

    def test_create_config_file(self, tmp_path):
        # configuration._override_config_file =
        with set_temp_config(tmp_path):
            assert dc.create_config_file("dysh", rootname="dysh-test")
            path = dc.get_config_dir_path(rootname="dysh-test")
        print(tmp_path)
        assert (tmp_path / "dysh-test/dysh.cfg").is_file()
        # Clean up.
        # (path / "dysh.cfg").unlink()
        # path.rmdir()
