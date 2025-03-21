"""Tests for dysh.config.core"""

import dysh.config as dc


class TestConfig:

    def test_get_config_dir_path(self):
        path = dc.get_config_dir_path(rootname="dysh-test")
        assert path.exists()

    def test_create_config_file(self):
        assert dc.create_config_file("dysh", rootname="dysh-test")
        path = dc.get_config_dir_path(rootname="dysh-test")
        print(path)
        assert (path / "dysh.cfg").is_file()
        # Clean up.
        (path / "dysh.cfg").unlink()
        path.rmdir()
