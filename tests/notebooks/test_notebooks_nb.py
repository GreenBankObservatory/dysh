
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

from dysh import util

root_dir = util.get_project_root()

def test_positionswitch():
    """
    """

    notebook = f"{root_dir}/notebooks/examples/positionswitch-test.ipynb"
    filename = f"{root_dir}/testdata/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"

    with open(notebook) as f:
        nb = nbformat.read(f, as_version=4)

    has_filename = False

    for cell in nb["cells"]:
        try:
            if "wget" in cell["metadata"]["tags"]:
                has_filename = True
                cell["source"] = f"filename = '{filename}'"
        except KeyError:
            continue

    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    exec_nb, resources = ep.preprocess(nb, {'metadata': {'path': '.'}})

    print_test = False

    for cell in exec_nb["cells"]:
        try:
            if "print_tsys" in cell["metadata"]["tags"]:
                assert cell["outputs"][0]["text"] == 'T_sys = 17.24\n'
                print_test = True
        except KeyError:
            continue

    assert print_test
