from astropy.utils.data import download_file as astropy_download_file


def download_file(*args, **kwargs):
    kwargs.setdefault("pkgname", "dysh")
    kwargs.setdefault("cache", True)
    return astropy_download_file(*args, **kwargs)
