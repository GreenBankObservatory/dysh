import os
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

# Collect notebook files at module level
NOTEBOOK_DIR = Path("notebooks/examples/")
NOTEBOOK_FILES = list(NOTEBOOK_DIR.glob("*.ipynb"))
ALLOWED_ERROR_NAMES = [
    "requests.exceptions.ReadTimeout",
    "requests.exceptions.HTTPError",
    "requests.exceptions.ConnectTimeout",
]

# Notebooks that require internet access (e.g. query external services like Splatalogue)
INTERNET_REQUIRED_NOTEBOOKS = {"line_search.ipynb"}


def _notebook_params():
    params = []
    for f in NOTEBOOK_FILES:
        marks = [pytest.mark.notebooks, pytest.mark.slow]
        if f.name in INTERNET_REQUIRED_NOTEBOOKS:
            marks.append(pytest.mark.requires_internet)
        params.append(pytest.param(f, marks=marks))
    return params


def check_notebook_execution(notebook_file):
    """Execute a given notebook and check for errors"""

    # Get the current testing directory
    cwd = Path().cwd()
    # Get the directory of the notebook
    notebook_dir = Path(notebook_file).parent.resolve()

    # Open the notebook
    with open(notebook_file, encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    # Try to execute the notebook
    try:
        # Change the current working directory to the notebook's directory
        os.chdir(notebook_dir)
        client = NotebookClient(nb, timeout=600, allow_error_names=ALLOWED_ERROR_NAMES)
        client.execute()
    finally:
        # Always restore directory
        os.chdir(cwd)


@pytest.mark.parametrize("notebook_file", _notebook_params(), ids=lambda x: x.name)
def test_notebook_execution(notebook_file):
    """Test that each notebook runs without errors"""
    check_notebook_execution(notebook_file)


@pytest.mark.notebooks
@pytest.mark.parametrize("notebook_file", NOTEBOOK_FILES, ids=lambda x: x.name)
def test_single_top_level_header(notebook_file):
    """Test that each notebook has exactly 1 top-level header"""

    with open(notebook_file, encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    top_level_headers = 0

    # Search each markdown cell for top-level headers
    for cell in nb.cells:
        if cell.cell_type == "markdown":
            source = cell.source
            lines = source.split("\n")
            md_code = False
            for line in lines:
                if line.startswith("```") and not md_code:
                    md_code = True
                elif line.startswith("```") and md_code:
                    md_code = False
                if line.startswith("# ") and not md_code:
                    top_level_headers += 1

    assert top_level_headers == 1, f"Notebook {notebook_file} has {top_level_headers} top-level headers."
