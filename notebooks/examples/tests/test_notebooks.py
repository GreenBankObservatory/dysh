import nbformat
from nbclient import NotebookClient
from glob import glob
import os
from pathlib import Path
import pytest

def check_notebook_execution(notebook_file):
    """Execute a given notebook and check for errors"""

    # Get the current testing directory
    cwd = Path().cwd()
    # Get the directory of the notebook
    notebook_dir = Path(notebook_file).parent.resolve()

    # Open the notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
    
    # Try to execute the notebook
    try:
        # Change the current working directory to the notebook's directory
        os.chdir(notebook_dir)
        client = NotebookClient(nb, timeout=600)
        client.execute()
        # Return to original working directory
        os.chdir(cwd)
    except Exception as e:
        # Return to original working directory
        os.chdir(cwd)
        pytest.fail(f"Error executing notebook {notebook_file}: {e}")

def test_notebooks_execution():
    """Test that each notebook runs without errors"""
    # List of notebook files to test
    notebook_files = glob("notebooks/examples/*.ipynb")
    
    for notebook_file in notebook_files:
        check_notebook_execution(notebook_file)

def test_single_top_level_header():
    """Test that each notebook has exactly 1 top-level header"""
    # List of notebook files to test
    notebook_files = glob("notebooks/examples/*.ipynb")
    
    # Open each notebook
    for notebook_file in notebook_files:
        with open(notebook_file, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
        
        top_level_headers = 0
        
        # Search each markdown cell for top-level headers
        for cell in nb.cells:
            if cell.cell_type == 'markdown':
                source = cell.source
                lines = source.split('\n')
                for line in lines:
                    if line.startswith('# '):
                        top_level_headers += 1

        assert top_level_headers == 1, f"Notebook {notebook_file} has {top_level_headers} top-level headers."

if __name__ == "__main__":
    test_notebooks_execution()
