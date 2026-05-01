import functools

import pytest
import requests

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
