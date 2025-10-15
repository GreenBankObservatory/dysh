"""Tests for dysh.config.core"""

from astropy.config import set_temp_config

import dysh.config as dc


class TestConfig:
    def test_get_config_dir_path(self):
        path = dc.get_config_dir_path(rootname="dysh-test")
        assert path.exists()

    def test_create_config_file(self, tmp_path):
        with set_temp_config(tmp_path):
            assert dc.create_config_file("dysh", rootname="dysh-test")
            path = dc.get_config_dir_path(rootname="dysh-test")  # noqa: F841
        assert (tmp_path / "dysh-test/dysh.cfg").is_file()

    def test_create_config_file_overwrite(self, tmp_path):
        conf_str = """[fits]\n## Maximum number of rows to be displayed by summary\nsummary_max_rows = 2"""
        with set_temp_config(tmp_path):
            assert dc.create_config_file("dysh", rootname="dysh-test")
            path = dc.get_config_dir_path(rootname="dysh-test")  # noqa: F841
            with open(tmp_path / "dysh-test/dysh.cfg", "w") as log:
                log.write(conf_str)
            with open(tmp_path / "dysh-test/dysh.cfg") as log:
                lines = log.read()
                assert lines == conf_str
            dc.create_config_file("dysh", rootname="dysh-test")
            with open(tmp_path / "dysh-test/dysh.cfg") as log:
                lines = log.read()
                assert lines == conf_str
