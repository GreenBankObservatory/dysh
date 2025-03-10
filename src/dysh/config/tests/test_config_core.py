"""Tests for dysh.config.core"""

import dysh.config as dc


class TestConfig:

    def test_get_config_dir_path(self):
        path = dc.get_config_dir_path(rootname="test-dysh")
        assert path.exists()

    def test_create_config_file(self):
        assert dc.create_config_file("dysh", rootname="test-dysh")
        # Clean up.
        path = dc.get_config_dir_path(rootname="test-dysh")
        (path / "dysh.cfg").unlink()
        path.rmdir()
