import wget


def test_FITS_wget():
    url = "https://www.gb.nrao.edu/dysh/examples/nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.A.fits"
    try:
        fname = wget.download(url)
    except Exception as exc:
        assert False, f"test_FITS_wget raised the following exception: {exc}"

def test_flag_wget():
    url = "https://www.gb.nrao.edu/dysh/examples/nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.A.flag"
    try:
        fname = wget.download(url)
    except Exception as exc:
        assert False, f"test_flag_wget raised the following exception: {exc}"