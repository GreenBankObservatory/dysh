import wget
import os

def test_flag_wget():
    url = "https://www.gb.nrao.edu/dysh/example_data/nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.A.flag"
    try:
        fname = wget.download(url)
        os.remove(fname)
    except Exception as exc:
        assert False, f"test_flag_wget raised the following exception: {exc}"