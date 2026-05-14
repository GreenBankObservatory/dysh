"""
Regression tests for the history-logging decorators in :mod:`dysh.log`.

See https://github.com/GreenBankObservatory/dysh/issues/1108
"""

from dysh.spectra import Spectrum


def _smooth_history_entries(spec):
    return [h for h in spec._history if "smooth(" in h]


def test_log_call_to_history_positional_args_recorded():
    """Positional args (including strings) are recorded with quotes preserved."""
    s = Spectrum.fake_spectrum()
    s.smooth("gaussian", 5)
    entries = _smooth_history_entries(s)
    assert len(entries) == 1
    assert "smooth('gaussian',5,)" in entries[0]


def test_log_call_to_history_keyword_args_recorded():
    """Keyword args are recorded even when the function has no ``**kwargs``."""
    s = Spectrum.fake_spectrum()
    s.smooth(kernel="gaussian", width=5)
    entries = _smooth_history_entries(s)
    assert len(entries) == 1
    assert "smooth(kernel='gaussian',width=5,)" in entries[0]


def test_log_call_to_result_keyword_args_recorded():
    """The result's _history records keyword args (covers @log_call_to_result paths)."""
    s = Spectrum.fake_spectrum()
    q = s.smooth(kernel="gaussian", width=5)
    entries = _smooth_history_entries(q)
    assert any("smooth(kernel='gaussian',width=5,)" in h for h in entries)
