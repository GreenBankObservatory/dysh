import os
import nbformat
import nbclient


from . import nb_test
from dysh import util

root_dir = util.get_project_root()
test_dir = util.get_project_testdata()

# Hardcode the tag we will use to identify
# the cell that defines which file will be
# loaded.
FILENAME_TAG = "wget"


def test_positionswitch():
    """
    """

    # Which notebook we will test?
    notebook = os.path.join(root_dir, r"notebooks", r"examples", r"positionswitch-test.ipynb")
    # What file will we use during the test?
    filename = os.path.join(test_dir, r"TGBT21A_501_11", r"TGBT21A_501_11.raw.vegas.fits")

    with open(notebook) as f:
        nb = nbformat.read(f, as_version=4)

    nb_tester = nb_test.NotebookTester(nb)

    # Find the cell that defines the file to be loaded
    # and modify it to point to a file in the codebase.
    idx = nb_tester.get_cell_idx_from_tag(FILENAME_TAG)
    assert len(idx) == 1
    nb_tester.nb.cells[idx[0]]["source"] = f"filename = r'{filename}'"

    # Run all cells and keep the kernel client alive.
    nb_tester.execute(cleanup_kc=False)

    # Get the value of the variable `filename`
    # and check that it matches the `filename`
    # we defined for the test.
    ret_filename = nb_tester.get_ref("filename")
    assert ret_filename == filename

    # We can also do things like this to test values.
    tsys = nb_tester.get_ref("psscan[0].tsys.mean()")
    assert tsys == 17.240003306306875

    # We can also capture the cell output and use it
    # for testing.
    idx = nb_tester.get_cell_idx_from_tag("print_tsys")
    assert len(idx) == 1 # make sure there is only one cell with the tag.
    assert nb_tester.nb.cells[idx[0]]["outputs"][0]["text"] == 'T_sys = 17.24\n'

    # Cleanup.
    nb_tester._cleanup_kernel()
