"""Tests for LabGUI ipympl warning suppression (ipympl#488)."""

import warnings
from types import SimpleNamespace

import matplotlib.figure
import pytest


@pytest.fixture
def plotbase_stub():
    """Create a minimal object that satisfies LabGUI's plotbase interface."""
    stub = SimpleNamespace()
    stub.figure = matplotlib.figure.Figure()
    stub.figure.subplots(nrows=1, ncols=1)
    return stub


def _collect_traitlets_deprecations(caught):
    """Filter a list of recorded warnings down to traitlets DeprecationWarnings."""
    return [w for w in caught if issubclass(w.category, DeprecationWarning) and "traitlets" in str(w.filename)]


def test_ipympl_canvas_emits_traitlets_warning(plotbase_stub):
    """ipympl Canvas/FigureManager emit a traitlets DeprecationWarning (the upstream bug exists)."""
    from ipympl.backend_nbagg import Canvas, FigureManager

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        canvas = Canvas(plotbase_stub.figure)
        canvas.manager = FigureManager(canvas, 0)

    assert len(_collect_traitlets_deprecations(caught)) > 0, (
        "Expected at least one traitlets DeprecationWarning from ipympl"
    )


def test_labgui_suppresses_traitlets_warning(plotbase_stub):
    """LabGUI (dysh-lab entry point) suppresses the traitlets DeprecationWarning."""
    from dysh.plot.labgui import LabGUI

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        LabGUI(plotbase_stub)

    assert _collect_traitlets_deprecations(caught) == [], (
        f"LabGUI should suppress traitlets warnings, but got: "
        f"{[str(w.message) for w in _collect_traitlets_deprecations(caught)]}"
    )


def test_import_dysh_installs_traitlets_warning_filter():
    """Importing dysh should install a global filter for the traitlets warning.

    Run in a subprocess to avoid pytest's per-test warning filter management,
    which resets warnings.filters and would hide the dysh-installed filter.
    """
    import subprocess
    import sys

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import warnings, dysh; "
            "match = [f for f in warnings.filters "
            "if f[0] == 'ignore' and f[1] is not None "
            "and f[1].search('Passing unrecognized arguments to super')]; "
            "assert len(match) > 0, f'filter not found in {warnings.filters}'",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
