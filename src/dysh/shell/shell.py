import argparse
import sys
from pathlib import Path
from typing import List, Union

import IPython
from traitlets.config import Config

import dysh.fits
from dysh import __version__
from dysh.config import create_config_file
from dysh.fits.sdfitsload import SDFITSLoad
from dysh.log import init_logging

# TODO: Derive URLs from pyproject.toml?
BANNER = f"""--------------------------------------------------------------------------
                         Welcome to Dysh v{__version__}

    Example usage: https://dysh.readthedocs.io/
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
DEFAULT_COLORS = "Linux"


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
    parser.add_argument("-v", "--verbosity", help="Set logging verbosity", type=int, default=2, choices=[0, 1, 2, 3])
    parser.add_argument("--log", help="Specify log path", type=Path)
    parser.add_argument("-q", "--quiet", help="Silence DEBUG- and INFO-level logs to stderr", action="store_true")
    parser.add_argument("--version", help="Print version and exit", action="store_true")
    parser.add_argument("--skip-config", help="Skip creating a configuration file", action="store_true")
    return parser.parse_known_args()


def init_shell(
    *ipython_args,
    colors=DEFAULT_COLORS,
    profile: Union[str, Path] = "DEFAULT_PROFILE",
    sdfits_files=None,
    skip_config=False,
):
    c = Config()
    import numpy as np
    import pandas as pd
    from astropy.io import fits
    from astropy.table import Table

    from dysh.fits.gbtfitsload import GBTFITSLoad, GBTOffline, GBTOnline
    from dysh.util.files import dysh_data

    user_ns = {
        "pd": pd,
        "np": np,
        "GBTFITSLoad": GBTFITSLoad,
        "GBTOnline": GBTOnline,
        "GBTOffline": GBTOffline,
        "dysh_data": dysh_data,
        "Table": Table,
        "fits": fits,
    }

    c.BaseIPythonApplication.profile = profile
    c.InteractiveShell.colors = colors
    c.InteractiveShell.banner2 = BANNER.format(
        user_ns_str="\n".join(f"{' ' * 8}{k} (from {v.__name__})" for k, v in user_ns.items())
    )
    if sdfits_files:
        user_ns["sdfits_files"] = sdfits_files
    if not skip_config:
        create_config_file("dysh", rootname="dysh")
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
    if args.version:
        print(__version__)
        sys.exit(0)
    init_logging(verbosity=args.verbosity, path=args.log, quiet=args.quiet)
    sdfits_files = open_sdfits_files(args.paths, args.fits_loader)
    init_shell(
        *remaining_args,
        colors=args.colors,
        profile=args.profile,
        sdfits_files=sdfits_files,
        skip_config=args.skip_config,
    )


if __name__ == "__main__":
    main()
