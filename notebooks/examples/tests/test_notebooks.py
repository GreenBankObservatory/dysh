import os
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

# Collect notebook files at module level
NOTEBOOK_DIR = Path("notebooks/examples/")
NOTEBOOK_FILES = list(NOTEBOOK_DIR.glob("*.ipynb"))


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
        client = NotebookClient(nb, timeout=600)
        client.execute()
    finally:
        # Always restore directory
        os.chdir(cwd)


@pytest.mark.notebooks
@pytest.mark.slow
@pytest.mark.parametrize("notebook_file", NOTEBOOK_FILES, ids=lambda x: x.name)
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
            for line in lines:
                if line.startswith("# "):
                    top_level_headers += 1

    assert top_level_headers == 1, f"Notebook {notebook_file} has {top_level_headers} top-level headers."
