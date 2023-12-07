from glob import glob

from docstr_coverage import analyze


def test_docstring_coverage():
    """
    Checks that the percentage of Python functions with docstrings
    is >= 0%
    """
    pyfiles = glob("./**/*.py", recursive=True)
    coverage_report = analyze(pyfiles)
    coverage = coverage_report.count_aggregate().coverage()
    assert coverage >= 0
