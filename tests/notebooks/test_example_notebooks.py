
from testbook import testbook

from dysh import util

root_dir = util.get_project_root()

@testbook(f'{root_dir}/notebooks/examples/positionswitch.ipynb', execute=False)
def test_notebook(tb):
    """
    """

    tb.inject(
    """
    filename = "/home/dysh/example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits"
    """,
    after=2,
    run=True
            )
    tb.execute_cell(["import"])
    tb.execute_cell(["load"])
    tb.execute_cell(["summary"])
    tb.execute_cell(["test"])
    tb.execute_cell(["timeaverage"])

    #ta = tb.ref("ta")
    #print(ta[0].meta["TSYS"])
