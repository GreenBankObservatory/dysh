import functools

import pytest
import requests
from astropy.table import Table

import dysh.line.search as _search_module

SPLATALOGUE_NETWORK_EXCEPTIONS = (
    requests.exceptions.HTTPError,
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
)


def skip_if_splatalogue_down(func):
    """Decorator: convert Splatalogue network errors to pytest.skip."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except SPLATALOGUE_NETWORK_EXCEPTIONS as e:
            pytest.skip(f"Splatalogue unavailable: {e}")

    return wrapper


@pytest.fixture
def mock_splatalogue_query(monkeypatch):
    """Patch Splatalogue.query_lines to return a predictable two-row table.

    The table mimics the raw Splatalogue response: orderedfreq in MHz as plain
    floats and intintensity as strings (as the real service returns them).
    """

    def _make_table(*args, **kwargs):
        return Table(
            {
                "orderedfreq": [115271.202, 230538.000],  # CO J=1-0, J=2-1 in MHz
                "name": ["CO v=0 1-0", "CO v=0 2-1"],
                "intintensity": ["-4.3", "-3.7"],
            }
        )

    monkeypatch.setattr(_search_module.Splatalogue, "query_lines", _make_table)
