import pytest

from dysh import util


@pytest.fixture(scope="module")
def root_dir():
    return util.get_project_root()


@pytest.fixture(scope="module")
def data_dir(root_dir):
    return util.get_project_testdata()
