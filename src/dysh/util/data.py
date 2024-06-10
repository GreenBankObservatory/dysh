""" """

import os
import tarfile
from pathlib import Path

import wget


def get_example_data(example, destination=None, rootname="dysh"):
    """ """

    dysh_url = "https://www.gb.nrao.edu/dysh/example_data/"

    examples = {
        "positionswitch": {
            "url": f"{dysh_url}/positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits",
            "action": None,
        },
        "frequencyswitch": {
            "url": f"{dysh_url}/frequencyswitch/data/TREG_050627/TREG_050627.raw.acs/TREG_050627.raw.acs.fits",
            "action": None,
        },
        "subbeamnod": {
            "url": f"{dysh_url}/subbeamnod/data/AGBT13A_124_06/AGBT13A_124_06.raw.acs/AGBT13A_124_06.raw.acs.fits",
            "action": None,
        },
        "multiplefiles": {"url": f"{dysh_url}/multifits/data/TRCO_230413_Ka.raw.vegas.tgz", "action": unpack},
    }

    if example not in examples.keys():
        raise KeyError(f"{example} is not available. The following are: {', '.join(examples.keys())}")

    # Handle destination dir.
    if destination is None:
        xch = os.environ.get("XDG_CACHE_HOME")
        # Use `XDG_CACHE_HOME` if it exists and is set.
        if xch is not None and os.path.exists(xch):
            xchpth = os.path.join(xch, rootname)
            # Check for dysh subdir.
            if os.path.exists(xchpth):
                dest = Path(xchpth)
            # Make it if not there already.
            else:
                dest = Path(xchpth)
                dest.mkdir()
        # Use the current directory if `XDG_CACHE_HOME` is not set or does not exist.
        else:
            dest = Path(".")
    # User provided destination.
    else:
        dest = Path(destination)

    if os.path.exists(dest):
        output = dest / examples[example]["url"].split("/")[-1]
        # Check if file exsists.
        if output.is_file():
            filename = output
        else:
            print(f"Will download data from: {examples[example]['url']}")
            filename = Path(wget.download(examples[example]["url"], out=str(output)))

    if examples[example]["action"] is not None:
        filename = examples[example]["action"](filename, filename.parent)

    print("")
    print(f"Data is available in : {filename}")

    return filename


def unpack(filename, dest):
    """ """

    tar = tarfile.open(filename)
    tar.extractall(path=dest)

    return dest.joinpath(filename.stem)
