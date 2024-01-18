import argparse
from pathlib import Path
from typing import List, Union

import IPython
from traitlets.config import Config

import dysh.fits
from dysh import __version__
from dysh.fits.sdfitsload import SDFITSLoad

# TODO: Derive URLs from pyproject.toml?
BANNER = f"""--------------------------------------------------------------------------
                         Welcome to Dysh v{__version__}

    Example usage: https://dysh.readthedocs.io/en/latest/example.html
    Bug reports:    https://github.com/GreenBankObservatory/dysh/issues

    For help with a Dysh routine from the command line,
    use the builtin 'help'. e.g.:

        help(GBTFITSLoad)

    The following items have been made globally available for convenience:
{{user_ns_str}}

   To suppress this banner, start with --no-banner flag
--------------------------------------------------------------------------
"""

DEFAULT_PROFILE = "dysh"
DEFAULT_COLORS = "LightBG"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Dysh interactive shell.\n\n All CLI arguments other than those defined below are passed through "
            "to ipython; see $ ipython --help for more details"
        )
    )
    parser.add_argument("paths", help="FITS file paths to load initially", nargs="*", type=Path)
    parser.add_argument("-p", "--profile", help="The IPython profile to use", default=DEFAULT_PROFILE)
    parser.add_argument("-L", "--fits-loader", help="The SDFITS loader class name to use", default="GBTFITSLoad")
    parser.add_argument(
        "--colors",
        help="Set the color scheme",
        choices=["NoColor", "Neutral", "Linux", "LightBG"],
        default=DEFAULT_COLORS,
    )
    return parser.parse_known_args()


def init_shell(*ipython_args, colors=DEFAULT_COLORS, profile: Union[str, Path] = "DEFAULT_PROFILE", sdfits_files=None):
    c = Config()
    import numpy as np
    import pandas as pd
    from astropy.io import fits
    from astropy.table import Table

    from dysh.fits.gbtfitsload import GBTFITSLoad

    user_ns = {"pd": pd, "np": np, "GBTFITSLoad": GBTFITSLoad, "Table": Table, "fits": fits}
    if sdfits_files:
        user_ns["sdfits_files"] = sdfits_files

    c.BaseIPythonApplication.profile = profile
    c.InteractiveShell.colors = colors
    c.InteractiveShell.banner2 = BANNER.format(
        user_ns_str="\n".join(f"{' '*8}{k} (from {v.__name__})" for k, v in user_ns.items())
    )
    IPython.start_ipython(ipython_args, config=c, user_ns=user_ns)


def get_fits_loader_class(loader_class_name: str):
    try:
        return getattr(dysh.fits, loader_class_name)
    except AttributeError as error:
        raise NotImplementedError(f"No known SDFITS Loader {loader_class_name!r}") from error


def open_sdfits_files(paths: List[Path], loader_class_name="GBTFITSLoad") -> List[SDFITSLoad]:
    loader_class = get_fits_loader_class(loader_class_name)
    return [loader_class(path) for path in paths]


def main():
    args, remaining_args = parse_args()
    sdfits_files = open_sdfits_files(args.paths, args.fits_loader)
    init_shell(*remaining_args, colors=args.colors, profile=args.profile, sdfits_files=sdfits_files)


if __name__ == "__main__":
    main()
