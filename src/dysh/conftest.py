"""
pytest configuration file for src/dysh.
"""


def pytest_collection_modifyitems(items):
    """
    Reorder the tests, so that shell test are last.
    Running shell tests before test that plot, breaks something.
    This was the solution I could find.

    See:
    https://docs.pytest.org/en/latest/how-to/writing_hook_functions.html#writinghooks
    and
    https://docs.pytest.org/en/latest/how-to/writing_plugins.html#writing-plugins
    for more details.
    """
    reordered = items[:]
    for item in items:
        if "shell" in item.name:
            reordered.remove(item)
            reordered.append(item)

    items[:] = reordered
