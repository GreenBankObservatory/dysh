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
from dysh.util.download import from_url

from ..util import minimum_string_match

_debug = False
# _debug = True

# note the examples in https://gbtdocs.readthedocs.io/en/latest/how-tos/data_reduction/gbtidl.html
# @todo   convert everything to use Path()
#         Path() cannot be used on input.... input needs to be a string

# fmt:off

# $DYSH/testdata
valid_dysh_test = {
    "test1"      : "AGBT05B_047_01/AGBT05B_047_01.raw.acs/",   # same as example='test1'
    "getps"      : "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits",
    "getfs"      : "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits",
}


# http://www.gb.nrao.edu/dysh/example_data or /home/dysh/example_data or $DYSH_DATA/example_data
# @todo   see if we want the staff training datasets in here
valid_dysh_example = {
    "test1"      : "positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits",   # staff training PS      same as test='test1'
    "getps"      : "onoff-L/data/TGBT21A_501_11.raw.vegas.fits",
                   #    positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/"
    "getpslarge" : "onoff-L/data/TGBT21A_501_11.raw.vegas/",
    "getfs"      : "fs-L/data/AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits",
                   #    frequencyswitch/data/TREG_050627/TREG_050627.raw.acs/"    # staff training FS
    "subbeamnod" : "subbeamnod-Ka/data/TRCO_230413_Ka.raw.vegas/TRCO_230413_Ka.raw.vegas.A.fits",
                   #    subbeamnod/data/AGBT13A_124_06/AGBT13A_124_06.raw.acs/"   # staff training SBN
    "nod"        : "nod-KFPA/data/TGBT22A_503_02.raw.vegas",    # nodding example (scan 62,63)
                   #              TGBT22A_503_02.raw.vegas      # FS example in data_reduction (scan 64)
    "otf1"       : "mapping-L/data/TGBT17A_506_11.raw.vegas",   # OTF L-band example NGC6946
}


# /home/dysh/acceptance_testing or $DYSH_DATA/acceptance_testing
# in acceptance_testing/data
# AGBT05B_047_01  AGBT15B_244_07  AGBT18A_503_02  AGBT19A_473_41  TGBT18A_500_06
# AGBT13A_240_03  AGBT16B_392_01  AGBT18B_014_02  AGBT19B_096_08  TGBT21A_501_10
# AGBT14B_480_06  AGBT17B_004_14  AGBT18B_354_03  AGBT20B_336_01  TREG_050627
# AGBT15B_228_08  AGBT17B_319_06  AGBT19A_080_01  AGBT22A_325_15  TSCAL_19Nov2015
valid_dysh_accept = {
    "nod1"       : "AGBT22A_325_15/AGBT22A_325_15.raw.vegas",
    "nod2"       : "TREG_050627/TREG_050627.raw.acs/TREG_050627.raw.acs.fits",               # deprecated
    "nod3"       : "AGBT15B_244_07/AGBT15B_244_07.raw.vegas",
    "nod4"       : "TGBT18A_500_06/TGBT18A_500_06.raw.vegas",
    "nod5"       : "TSCAL_19Nov2015/TSCAL_19Nov2015.raw.acs/TSCAL_19Nov2015.raw.acs.fits",   # deprecated
    "nod6"       : "AGBT17B_319_06/AGBT17B_319_06.raw.vegas",
    "nod7"       : "TGBT21A_501_10/TGBT21A_501_10.raw.vegas",
    "nod8"       : "AGBT19A_340_07/AGBT19A_340_07.raw.vegas",
    "nod9"       : "AGBT12A_076_05/AGBT12A_076_05.raw.acs",
    "multismallsmall" : "AGBT20B_336_01/AGBT20B_336_01.raw.vegas",  # multiple small FITS files (54M each), small flags files (7 lines)
    "multihugesmall"  : "AGBT14B_480_06/AGBT14B_480_06.raw.vegas",  # multiple huge FITS files (3.5GM each), small flags files (6 lines)
    "multismallbig" : "AGBT23A_432_03/AGBT23A_432_03.raw.vegas", # multiple small FITS files (64M each), large flag files (20 lines)
    "multibighuge"  : "AGBT17B_319_06/AGBT17B_319_06.raw.vegas",  # multiple large FITS files (733M each), huge flag files (102 lines)

}

# fmt: on


def dysh_data(sdfits=None, test=None, example=None, accept=None, dysh_data=None, verbose=False):
    r"""
    Resolves the filename within the GBO dysh data system without the need for an absolute path.

    Currently configured to work at GBO, where for example /home/sdfits exists. For other sites users
    need to configure a $DYSH_DATA directory, properly populated with (symlinks to) project and test data,
    as described below. Optionally, an explicit dysh_data= can be given, which overrides any possible $DYSH_DATA
    environment (or configuration) that may exist.

    Only one of the keywords sdfits=, test=, example=, accept= can be given to probe for data. They are
    processed in that order, whichever comes first.


    Locations of various dysh_data directory roots:  ($DYSH is the repo root for developers)
    -----------------------------------------------
    keyword       location at GBO                      $DYSH_DATA root
    -------       ---------------                      ---------------
    sdfits=       /home/sdfits                         $DYSH_DATA/sdfits
    test=         $DYSH/testdata                       $DYSH_DATA/testdata
    example=      /home/dysh/example_data              $DYSH_DATA/example_data
    accept=       /home/dysh/acceptance_testing        $DYSH_DATA/acceptance_testing

    Note: test= resolves to the same filename as the util.get_project_testdata() function


    Examples of use including mnemonics or full paths:
    --------------------------------------------------
    fn = dysh_data(test='getps')
    fn = dysh_data(example='getfs')
    fn = dysh_data(example='onoff-L/data/TGBT21A_501_11.raw.vegas')
    fn = dysh_data('AGBT21B_024_54')         ->  /home/sdfits/AGBT21B_024_54
                                        or:  -   /lma1/teuben/GBT-EDGE/rawdata/AGBT21B_024_54 with $DYSH_DATA


    Notes:

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

    def sdfits_offline(fn):
        """fn is an sdfits= filename that was shown to exist
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
    if verbose:
        print("DYSH_DATA:", dysh_data)

    # Override the default _url.
    if "DYSH_DATA_URL" in os.environ:
        _url = os.environ["DYSH_DATA_URL"]

    # 2. Process whichever one of 'sdfits=', 'test=', 'example=', and  'accept=' is present (in that order)

    # sdfits:   the main place where GBO data reside

    if sdfits != None:  # noqa: E711
        if sdfits == "?":
            if dysh_data == None:  # noqa: E711
                dd = Path("/home/sdfits")
            else:
                dd = dysh_data / "sdfits"
            # @todo figure out listing of file OS agnostic
            cmd = "ls %s" % dd
            print("# dysh_data::sdfits")
            print("# contents of", dd)
            print("# -----------------")
            os.system(cmd)
            return None
        if dysh_data != None:  # noqa: E711
            fn = dysh_data / Path("sdfits") / sdfits  # normally user is using a private sdfits
            if fn.exists():
                return sdfits_offline(fn)
        fn = Path("/home/sdfits/") / sdfits  # expected at GBO
        if fn.exists():
            return sdfits_offline(fn)
        # print(f"could not handle sdfits={sdfits} yet")
        return None

    # test:   this should also be allowed to use util.get_project_testdata() as well

    if test != None:  # noqa: E711
        if test == "?":
            print("# dysh_data::test")
            print("# ---------------")
            for k in valid_dysh_test.keys():
                print(k, valid_dysh_test[k])
            return None
        my_test = minimum_string_match(test, list(valid_dysh_test.keys()))
        if my_test != None:  # noqa: E711
            my_test = valid_dysh_test[my_test]
        else:
            my_test = test
        #
        if dysh_data != None:  # noqa: E711
            fn = dysh_data / "testdata" / my_test
            if not fn.exists():
                fn = util.get_project_testdata() / my_test
        else:
            fn = util.get_project_testdata() / my_test
        if verbose:
            print("final:", fn)
        if fn.exists():  # @todo this catches files and directories
            return fn
        print("Could not find", fn)
        return None

    # example:  these can also obtain data via from_url (or perhaps astropy caching???)

    if example != None:  # noqa: E711
        if example == "?":
            print("# dysh_data::example")
            print("# ------------------")
            for k in valid_dysh_example.keys():
                print(k, valid_dysh_example[k])
            return None
        my_example = minimum_string_match(example, list(valid_dysh_example.keys()))
        if my_example != None:  # noqa: E711
            my_example = valid_dysh_example[my_example]
        else:
            my_example = example
        if dysh_data != None:  # noqa: E711
            fn = dysh_data / "example_data" / my_example
            if fn.exists():
                return fn
            print("Odd-1, did not find", fn)
        if dysh_data == None and os.path.exists(_example_data):  # noqa: E711
            fn = Path(_example_data) / my_example
            if fn.exists():
                return fn
            print("Odd-2, did not find", fn)
        # last resort, try getting it via from_url, but it will then be a local file in the current directory
        url = _url + "/example_data/" + my_example
        if verbose:
            print("url:", url)
        filename = url.split("/")[-1]
        if not os.path.exists(filename):
            print(f"Downloading {filename} from {url}")
            try:
                filename = from_url(url)
                print(f"\nRetrieved {filename}")
            except:  # noqa: E722
                print(f"\nFailing to retrieve example {filename} ")
                return None
        else:
            print(f"{filename} already downloaded")
        return Path(filename)

    # accept:   acceptance_testing/data - from_url not recommended (does not work on multifile fits)

    if accept != None:  # noqa: E711
        if accept == "?":
            print("# dysh_data::accept")
            print("# -----------------")
            for k in valid_dysh_accept.keys():
                print(k, valid_dysh_accept[k])
            return None
        my_accept = minimum_string_match(accept, list(valid_dysh_accept.keys()))
        if my_accept != None:  # noqa: E711
            my_accept = valid_dysh_accept[my_accept]
        else:
            my_accept = accept
        if dysh_data != None:  # noqa: E711
            fn = dysh_data / "acceptance_testing/data" / my_accept
            if fn.exists():
                return fn
            print("Odd-1, did not find", fn)
        if dysh_data == None and os.path.exists(_accept_data):  # noqa: E711
            fn = Path(_accept_data) / my_accept
            if fn.exists():
                return fn
            print("Odd-2, did not find", fn)
        # last resort, try getting it via from_url, but it will then be a local file in the current directory
        url = _url + "/acceptance_testing/data/" + my_accept
        if verbose:
            print("url:", url)
        filename = url.split("/")[-1]
        if not os.path.exists(filename):
            print(f"Downloading {filename} from {url}")
            try:
                filename = from_url(url)
                print(f"\nRetrieved {filename}")
            except:  # noqa: E722
                print(f"\nFailing to retrieve accept {filename}")
                return None
        else:
            print(f"{filename} already downloaded")
        return Path(filename)

    print("You have not given one of:   sdfits=, test=, example=, accept=")
    print("or use =? as argument to get a list of valid shortcuts")
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
        if _debug:
            print("# FNAME:", fname)

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
                    print("# Warning: %s not in the environment" % p)
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
                print("# Warning: directory %s does not exist" % p)
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
    p.add_argument("-m", "--maxfiles", type=int, default=None, help="Maximum number of files to return [Default: all]")
    p.add_argument("-c", "--count", action="store_true", help="add counter to filenames?")
    p.add_argument("-w", "--wildcard", action="store_true", help="fully wildcard the filename embedded")
    p.add_argument("-r", "--recursive", action="store_true", help="resursive?")
    p.add_argument("-p", "--path", type=str, default=None, help="optional (colon separated) path(s)")
    p.add_argument("filename", nargs="+", help="Filename(s) to search for")

    args = p.parse_args()
    if _debug:
        print("#", args)

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
