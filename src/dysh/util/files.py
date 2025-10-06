#! /usr/bin/env python
#
#  Tools to operate on files
#     dysh_data = simple frontend to grab common dysh filenames
#     fdr       = Recursive data (file) finder
#
#  command line usage:
#
#      fdr [-r] [-m MAXFILES] file [path1 path2 ...]
#


import glob
import os
from pathlib import Path

import dysh.util as util
from dysh.log import logger
from dysh.util.download import from_url

from ..util import minimum_string_match

# the GBTIDL examples from https://gbtdocs.readthedocs.io/en/latest/how-tos/data_reduction/gbtidl.html
#         getps:        "data/ngc2415.fits"                NGC2415    example=getps2    (TGBT21A_501_11)
#         getfs:        "data/TGBT22A_503_02.raw.vegas"    W3_1       test=    example=
#         getsigref:    "data/TGBT22A_503_02.raw.vegas"    W3_1
#         getps:        "data/AGBT17A_404_01.raw.vegas"    A123606    test=
#
# @todo   convert everything to use Path()
#         Path() cannot be used on input.... input needs to be a string

# fmt:off

# $DYSH/testdata      @ todo   normalize names with the example= cases
# ~300 MB
valid_dysh_test = {
    "getps"      : "AGBT05B_047_01/AGBT05B_047_01.raw.acs/",                                  # NGC5291   old test1
    "getps2"     : "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits",                            # NGC2415   one int, 540k
    "getfs"      : "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits", # W3OH
    "subbeamnod" : "TRCO_230413_Ka",

}


# http://www.gb.nrao.edu/dysh/example_data or /home/dysh/example_data or $DYSH_DATA/example_data
# @todo   see if we want the staff training datasets in here
# ~410 GB
valid_dysh_example = {
    "getps"      : "positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits", #  NGC5291  old test1
                   # Used in a lot of example notebooks:
                   #   example/dataIO
                   #   example/metadata_management
                   #   example/positionswitch
                   #   example/smoothing
                   #   example/velocity frames
    "getps0"     : "positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs",                             #  NGC5291  old test1
    "getps2"     : "onoff-L/data/TGBT21A_501_11.raw.vegas.fits",    #  NGC2415   - old getps
    "getpslarge" : "onoff-L/data/TGBT21A_501_11.raw.vegas/",        #  NGC2415, NGC2782 etc. - total 15GB
    "getfs"      : "fs-L/data/AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits",
    "getfs2"     : "frequencyswitch/data/TREG_050627/TREG_050627.raw.acs/",    #  W3OH    # staff training FS
    "subbeamnod" : "subbeamnod/data/AGBT13A_124_06/AGBT13A_124_06.raw.acs/",   #  vIIzw31      example/subbeamnod    staff training SBN -- no signal?
    "subbeamnod2": "subbeamnod-Ka/data/TRCO_230413_Ka.raw.vegas/TRCO_230413_Ka.raw.vegas.A.fits",
    "nod"        : "nod-KFPA/data/TGBT22A_503_02.raw.vegas/",       # W3_1      example/nodding  (scan 62,63)
                   #              TGBT22A_503_02.raw.vegas          # FS example in data_reduction (scan 64)
    "align"      : "mixed-fs-ps/data/TGBT24B_613_04.raw.vegas.trim.fits",  #   MESSIER32  example/align_spectra
    "flagging"   : "rfi-L/data/AGBT17A_404_01.tar.gz",                     # tar.gz not yet supported?     A123606  example/flagging
    "survey"     : "hi-survey/data/AGBT04A_008_02.raw.acs/AGBT04A_008_02.raw.acs.fits",   # example/hi-survey
    "otf1"       : "mapping-L/data/TGBT17A_506_11.raw.vegas/",      # OTF L-band NGC6946
    "otf2"       : "AGBT21B_024_01",                                # OTF Argus  NGC0001 (EDGE)
    "otf3"       : "AGBT21B_024_20",                                # OTF Argus  NGC5954 (EDGE)
    "otf4"       : "mapping-Argus/data/TGBT22A_603_05.raw.vegas/",  # OTF Argus  DR21
}


# /home/dysh/acceptance_testing or $DYSH_DATA/acceptance_testing
# in acceptance_testing/data
# AGBT05B_047_01  AGBT15B_244_07  AGBT18A_503_02  AGBT19A_473_41  TGBT18A_500_06
# AGBT13A_240_03  AGBT16B_392_01  AGBT18B_014_02  AGBT19B_096_08  TGBT21A_501_10
# AGBT14B_480_06  AGBT17B_004_14  AGBT18B_354_03  AGBT20B_336_01  TREG_050627
# AGBT1s5B_228_08  AGBT17B_319_06  AGBT19A_080_01  AGBT22A_325_15  TSCAL_19Nov2015
# ~ 33 GB
valid_dysh_accept = {
    "nod1"            : "AGBT22A_325_15/AGBT22A_325_15.raw.vegas",
    "nod2"            : "TREG_050627/TREG_050627.raw.acs/TREG_050627.raw.acs.fits",               # deprecated?   W3OH  example/frequencyswitch
    "nod3"            : "AGBT15B_244_07/AGBT15B_244_07.raw.vegas",                                # M82 examples/calseq
    "nod4"            : "TGBT18A_500_06/TGBT18A_500_06.raw.vegas",
    "nod5"            : "TSCAL_19Nov2015/TSCAL_19Nov2015.raw.acs/TSCAL_19Nov2015.raw.acs.fits",   # deprecated
    "nod6"            : "AGBT17B_319_06/AGBT17B_319_06.raw.vegas",
    "nod7"            : "TGBT21A_501_10/TGBT21A_501_10.raw.vegas",
    "nod8"            : "AGBT19A_340_07/AGBT19A_340_07.raw.vegas",
    "nod9"            : "AGBT12A_076_05/AGBT12A_076_05.raw.acs",
    "multismallsmall" : "AGBT20B_336_01/AGBT20B_336_01.raw.vegas",  # multiple small FITS files (54M each), small flags files (7 lines)
    "multihugesmall"  : "AGBT14B_480_06/AGBT14B_480_06.raw.vegas",  # multiple huge FITS files (3.5G each), small flags files (6 lines)
    "multismallbig"   : "AGBT23A_432_03/AGBT23A_432_03.raw.vegas",  # multiple small FITS files (64M each), large flag files (20 lines)
    "multibighuge"    : "AGBT17B_319_06/AGBT17B_319_06.raw.vegas",  # multiple large FITS files (733M each), huge flag files (102 lines)

}

# fmt: on


def dysh_data(sdfits=None, test=None, example=None, accept=None, dysh_data=None, gui=False):
    r"""Resolves the filename within the dysh data system without the need
    for an absolute path by passing mnemonics to any of four entry
    points (`sdfits`, `!test`, `example`, `accept`).

    Currently configured to work at GBO. For other sites users need to
    configure a $DYSH_DATA directory, properly populated with
    (symlinks to) directories as described below. Optionally, an
    explicit `dysh_data` can be given, which overrides any possible
    $DYSH_DATA environment (or configuration) that may exist.

    Only one of the keywords `sdfits`, `!test`, `example`, `accept` can be
    given to probe for data.

    As an exception, if the first argument (`sdfits`) has an absolute
    filename, it is passed unchecked.

    gui mode is experimental and may disappear or re-implemented at a
    later stage.


    The locations of various dysh_data directory roots are presented in the following
    Table, where $DYSH is the repo root for developers (this can be found using
    `dysh.util.get_project_root`).

    ========      =============================        =============================
    keyword       location at GBO                      $DYSH_DATA root
    ========      =============================        =============================
    sdfits=       /home/sdfits                         $DYSH_DATA/sdfits
    test=         $DYSH/testdata                       $DYSH_DATA/testdata
    example=      /home/dysh/example_data              $DYSH_DATA/example_data
    accept=       /home/dysh/acceptance_testing        $DYSH_DATA/acceptance_testing
    ========      =============================        =============================

    `!test` resolves to the same filename as the `util.get_project_testdata()` function
    but it is otherwise only available for developers (the testdata directory is not available
    if you `pip install dysh`).

    If present, the $SDFITS_DATA directory is honored instead of the default for `sdfits`
    and overrides the $DYSH_DATA directory.

    Examples
    --------

    Using mnemonics

    >>> fn = dysh_data(test='getps')
    >>> fn = dysh_data(example='getfs')

    Using full paths

    >>> fn = dysh_data(example='onoff-L/data/TGBT21A_501_11.raw.vegas')

    Using a project id

    >>> fn = dysh_data('AGBT21B_024_54')

    This will return `/home/sdfits/AGBT21B_024_54` at GBO, or `${DYSH_DATA}/sdfits/AGBT21B_024_54`
    if the $DYSH_DATA environment variable is set.

    Notes
    -----

    1) if $DYSH_DATA exist (and this is a new proposal), it will prepend
       that to the argument of get_dysh_data() and check for existence
       if $DYSH_DATA does not exist, but $SDFITS_DATA exists (a GBTIDL feature)
       it will use that
    2) if /home/dysh exists, it will prepend this and check for existence
       this will keep GBO people happy.  Offsite a symlink should also work.
    3) if none of those gave a valid name, it will fall back to making a URL
       by prepending http://www.gb.nrao.edu/dysh/ and using
       from_url for as long we want to support that.
       astropy caching is also an option
    4) directories (names not ending on .fits) cannot be downloaded using from_url
    5) configuration TBD

    """
    # fmt:off
    _url                = "http://www.gb.nrao.edu/dysh/"            # base of all things dysh
    _example_data       = "/home/dysh/public_html/example_data"     # GBO direct access
    _test_data          = "/home/dysh/public_html/test_data"        # not used ??
    _accept_data        = "/home/dysh/acceptance_testing/data"      # not in public_html ??
    # fmt:on

    if type(dysh_data) is str:
        dysh_data = Path(dysh_data)

    def sdfits_offline(fn):
        """fn is an sdfits= file or directory that was shown to exist
        If fn contains only one name
        See also GBTOffline()
        """
        if fn.is_file():
            return fn
        if not fn.is_dir():
            print(f"{fn} is not a file nor a directory, dunno how to proceed")
            return None
        # find all fits files one level deep
        ff = list(fn.glob("*/*.fits"))
        if len(ff) == 0:
            return fn
        # ensure there is only a single parent
        parents = []
        for f in ff:
            parents.append(f.parent)
        parents = list(set(parents))
        if len(parents) > 1:
            print(f"{fn} does not contain a single fits tree: {parents}")
            # @todo throw ?  or return the first one?

        return parents[0]

    def use_gui(my_dir):
        logger.debug(f"Using the GUI on {my_dir} is totally experimental.")
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename(initialdir=my_dir)
        # can currently only ask for files, use askdirectory() otherwise
        return file_path

    # 1.  find out if there is a dysh_data (or use $DYSH_DATA, or a .dyshrc config?)
    #     - if present, API dysh_data is used
    #     - if present, $DYSH_DATA is used
    #     - if present, python_env is used
    #     - if all of this fails, assume we're at GBO (all via /home/dysh)
    #     - if that still fails, look at current working directory
    #     - throw!?
    #     ? e.g. dysh_data('foo.fits') ->   sdfits='foo.fits'

    if dysh_data == None and "DYSH_DATA" in os.environ:  # noqa: E711
        dysh_data = Path(os.environ["DYSH_DATA"])
    logger.debug(f"DYSH_DATA: {dysh_data}")

    # 2. Process whichever one of 'sdfits=', 'test=', 'example=', and  'accept=' is present (in that order)

    # sdfits:   the main place where GBO data reside

    if sdfits is not None:
        if sdfits == "!":
            logger.warning("The GUI is experimental, it can only select a single fits file, no directories")
            return use_gui(dysh_data)

        if sdfits == "?" or sdfits == "*":
            if "SDFITS_DATA" in os.environ:
                dd = Path(os.environ["SDFITS_DATA"])
            elif dysh_data == None:  # noqa: E711
                dd = Path("/home/sdfits")
            else:
                dd = Path(dysh_data) / "sdfits"
            # @todo figure out listing of file OS agnostic
            cmd = f"ls {dd}"
            print("# dysh_data::sdfits")
            print("# contents of", dd)
            print("# -----------------")
            os.system(cmd)
            return None
        if dysh_data is not None:
            fn = dysh_data / Path("sdfits") / sdfits  # normally user is using a private sdfits
            if fn.exists():
                return sdfits_offline(fn)
        if "SDFITS_DATA" in os.environ:
            fn = Path(os.environ["SDFITS_DATA"]) / sdfits
            return sdfits_offline(fn)
        fn = Path("/home/sdfits/") / sdfits  # expected at GBO
        if fn.exists():
            return sdfits_offline(fn)
        # print(f"could not handle sdfits={sdfits} yet")
        return None

    # test:   this should also be allowed to use util.get_project_testdata() as well

    if test is not None:
        if test == "?":
            print("# dysh_data::test")
            print("# ---------------")
            for k in valid_dysh_test.keys():
                print(k, valid_dysh_test[k])
            return None
        my_test = minimum_string_match(test, list(valid_dysh_test.keys()))
        if my_test is not None:
            my_test = valid_dysh_test[my_test]
        else:
            my_test = test
        if dysh_data is not None:
            fn = dysh_data / "testdata" / my_test
            if not fn.exists():
                fn = util.get_project_testdata() / my_test
        else:
            fn = util.get_project_testdata() / my_test
        logger.debug(f"final: {fn}")
        if fn.exists():  # @todo this catches files and directories
            return fn
        print("Could not find", fn)
        return None

    # example:  these can also obtain data via from_url (or perhaps astropy caching???)

    if example is not None:
        if example == "?":
            print("# dysh_data::example")
            print("# ------------------")
            for k in valid_dysh_example.keys():
                print(k, valid_dysh_example[k])
            return None
        my_example = minimum_string_match(example, list(valid_dysh_example.keys()))
        if my_example is not None:
            my_example = valid_dysh_example[my_example]
        else:
            my_example = example
        if dysh_data is not None:
            fn = dysh_data / "example_data" / my_example
            if fn.exists():
                return fn
            print("Odd-1, did not find", fn)
        if dysh_data is None and os.path.exists(_example_data):
            fn = Path(_example_data) / my_example
            if fn.exists():
                return fn
            print("Odd-2, did not find", fn)
        # last resort, try getting it via from_url, but it will then be a local file in the current directory
        url = _url + "/example_data/" + my_example
        logger.info(f"url: {url}")
        filename = url.split("/")[-1]
        if not os.path.exists(filename):
            print(f"Downloading {filename} from {url}")
            try:
                filename = from_url(url)
                print(f"\nRetrieved {filename}")
            except Exception as e:
                print(f"\nFailing to retrieve example {filename} ")
                print(e)
                return None
        else:
            print(f"{filename} already downloaded")
        return Path(filename)

    # accept:   acceptance_testing/data - from_url not recommended (does not work on multifile fits)

    if accept is not None:
        if accept == "?":
            print("# dysh_data::accept")
            print("# -----------------")
            for k in valid_dysh_accept.keys():
                print(k, valid_dysh_accept[k])
            return None
        my_accept = minimum_string_match(accept, list(valid_dysh_accept.keys()))
        if my_accept is not None:
            my_accept = valid_dysh_accept[my_accept]
        else:
            my_accept = accept
        if dysh_data is not None:
            fn = dysh_data / "acceptance_testing/data" / my_accept
            if fn.exists():
                return fn
            print("Odd-1, did not find", fn)
        if dysh_data is None and os.path.exists(_accept_data):
            fn = Path(_accept_data) / my_accept
            if fn.exists():
                return fn
            print("Odd-2, did not find", fn)
        # last resort, try getting it via from_url, but it will then be a local file in the current directory
        url = _url + "/acceptance_testing/data/" + my_accept
        logger.debug(f"url: {url}")
        filename = url.split("/")[-1]
        if not os.path.exists(filename):
            print(f"Downloading {filename} from {url}")
            try:
                filename = from_url(url)
                print(f"\nRetrieved {filename}")
            except Exception as e:
                print(f"\nFailing to retrieve accept {filename}")
                print(e)
                return None
        else:
            print(f"{filename} already downloaded")
        return Path(filename)

    print("You have not given one of:   sdfits=, test=, example=, accept=")
    print("or use =? as argument to get a list of valid shortcuts")
    print(f"DYSH_DATA = {dysh_data}")
    return None


# def find_data_recursively(filename, path=None, recursive=False, wildcard=False, maxfiles=None):
def fdr(filename, path=None, recursive=False, wildcard=False, maxfiles=None):
    """
    Input:
        filename - can be wildcard too (but see the wildcard option)
        path - optional. can be : separated, can start with $ if envvar
        recursive - recursively search: Default: False
        wildcard - automatically wildcard the filename: Default not used
        maxfiles - maximum number of files to be returns. Default: All

    Returns:
        list of found filenames,  with maxfiles entries if applicable.
        Note list could be empty.
        Note if multiple paths are given,  maxfiles is applied to each sublist

    See also:
        astropy's getdata ???
        pdrptry.pdrutils.get_testdata()
        astropy.utils.data.get_pkg_data_filenames

    Examples:
        fdr('ngc1234.fits')    - this exact file!
        fdr('*.fits')          - all fits file in this directory
        fdr('ngc1234.fits','/tmp')  - this file in /tmp
        fdr('*.fits','/tmp')        - all fits files in /tmp
        fdr('ngc1234.fits','$DYSH_DATA_PATH')
        fdr('ngc1234.fits','$DYSH_DATA_PATH', True)
        fdr('ngc1234.fits','$DYSH_DATA_PATH:/data/gbt')
    """
    if os.path.exists(filename):
        return [filename]

    if path is None:
        if wildcard:
            fname = "*" + filename + "*"
        else:
            fname = filename
        if recursive:
            fname = "**/" + fname
        logger.debug("# FNAME:", fname)

        fn = glob.glob(fname, recursive=recursive)

        if maxfiles is None:
            retval = fn
        else:
            retval = fn[:maxfiles]
        retval.sort()
    else:
        cwd0 = os.getcwd()
        all = []
        for p in path.split(":"):
            if p[0] == "$":
                if p[1:] in os.environ:
                    p = os.environ[p[1:]]
                else:
                    print(f"# Warning: {p} not in the environment")
            if os.path.exists(p):
                os.chdir(p)
                if wildcard:
                    fname = "*" + filename + "*"
                else:
                    fname = filename
                if recursive:
                    fname = "**/" + fname
                fn = glob.glob(fname, recursive=recursive)
                fn.sort()
                if maxfiles is not None:
                    fn = fn[:maxfiles]
                all = all + fn
            else:
                print(f"# Warning: directory {p} does not exist")
            os.chdir(cwd0)
        retval = all
    return retval


def main_cli():
    import argparse

    my_help = """
    This script searches for files, optionally hierarchically,
    much like the Unix 'find' program.
    A difference is handling the --path directive, as multiple
    colon separated paths can be given, much like the $PATH
    environment variable in Unix.
    The path variable can also expand $-environment variables.

    """

    p = argparse.ArgumentParser(description=my_help, epilog="And so the search goes on....")
    p.add_argument(
        "-m",
        "--maxfiles",
        type=int,
        default=None,
        help="Maximum number of files to return [Default: all]",
    )
    p.add_argument("-c", "--count", action="store_true", help="add counter to filenames?")
    p.add_argument("-w", "--wildcard", action="store_true", help="fully wildcard the filename embedded")
    p.add_argument("-r", "--recursive", action="store_true", help="resursive?")
    p.add_argument("-p", "--path", type=str, default=None, help="optional (colon separated) path(s)")
    p.add_argument("filename", nargs="+", help="Filename(s) to search for")

    args = p.parse_args()
    logger.debug("#", args)

    filename = args.filename
    maxfiles = args.maxfiles
    recursive = args.recursive
    wildcard = args.wildcard
    path = args.path
    count = args.count

    r = []
    for f in filename:
        r = r + fdr(f, path, recursive, wildcard, maxfiles)
    # r.sort()

    if count:
        n = 1
        for f in r:
            print(n, f)
            n = n + 1
    else:
        for f in r:
            print(f)


if __name__ == "__main__":
    main_cli()
