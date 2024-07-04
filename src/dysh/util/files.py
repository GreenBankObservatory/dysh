#! /usr/bin/env python
#
#  Tools to operate on files
#     dysh_data = simple file grabber frontend
#     fdr       = Recursive data (file) finder
#
#  command line usage:
#
#      fdr [-r] [-m MAXFILES] file [path1 path2 ...]
#
#  Shell alternative:
#
#  find . -name \*$a\* -print


import os
import sys
import glob
import wget     # might get deprecated
from ..util import minimum_string_match
import dysh.util as util

_debug = True


# $DYSH/testdata 
valid_dysh_tests = {
    "getps" : "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits",
    "getfs" : "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits",
}


# http://www.gb.nrao.edu/dysh/example_data or /home/dysh/example_data
valid_dysh_examples = {
    "getps"      : "onoff-L/data/TGBT21A_501_11.raw.vegas.fits",
    "getfs"      : "fs-L/data/AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits",
    "subbeamnod" : "subbeamnod-Ka/data/TRCO_230413_Ka.raw.vegas/TRCO_230413_Ka.raw.vegas.A.fits",
}


def dysh_data(test=None,
              example=None,
              sdfits=None,
              url="http://www.gb.nrao.edu/dysh/",
              dysh_data = None):
    r"""

    Users that are not on the GBO system and need access to example_data could
    rsync the /home/dysh/example_data tree.  For example inside their $HOME/dysh_data/
    one would then set
             export DYSH_DATA=$HOME/dysh_data

    Locations of various dysh_data directory roots:
    -----------------------------------------------
    keyword       location                       method                           $DYSH_DATA root    
    test:         $DYSH/testdata                 util.get_project_testdata()      $DYSH_DATA/testdata
    example:      /home/dysh/example_data                                         $DYSH_DATA/example_data
    sdfits:       /home/sdfits                                                    $DYSH_DATA/sdfits


    Examples of use:
    ----------------
    fn = dysh_data(test='getps')
    fn = dysh_data(example='getfs')
    fn = dysh_data(sdfits='AGBT21B_024_54')         ->  /home/sdfits/AGBT21B_024_54
                                               or:  -   /lma1/teuben/GBT-EDGE/rawdata/AGBT21B_024_54


    Older Notes:

    0) if a local file AGBT20B_014_03.raw.vegas.A.fits already exist, 
       use it (thus supporting the old style)
    1) if $DYSH_DATA exist (and this is a new proposal), it will prepend 
       that to the argument of get_dysh_data() and check for existence
    2) if /home/dysh exists, it will prepend this and check for existence
       this will keep GBO people happy
    3) if none of those gave a valid name, it will fall back to making a URL 
       by prepending http://www.gb.nrao.edu/dysh/ and using
       wget for as long we want to support that. 
    """
    if dysh_data == None and 'DYSH_DATA' in os.environ:
        dysh_data = os.environ['DYSH_DATA']

    print("DYSH_DATA:", dysh_data)

    if test != None:
        my_test = minimum_string_match(test, list(valid_dysh_tests.keys()))
        fn = valid_dysh_tests[my_test]
        if dysh_data != None:
            fn = dysh_data + '/testdata/' + fn
            print("final-1:",fn)
            if os.path.exists(fn):    # @todo this catches files and directories
                return fn
            # weird, typo is need to check the code tree
        fn = util.get_project_testdata() / fn
        print('final-2:',fn)
        if os.path.exists(fn):    # @todo this catches files and directories        
            return fn
        print("Could not find file",fn)
        return None
        

    if example != None:
        my_example = minimum_string_match(example, list(valid_dysh_examples.keys()))
        fn = valid_dysh_examples[my_example]
        print('example:',fn)
        if dysh_data != None:
            fn1 = dysh_data + '/example_data/' + fn
            print("final-1",fn1)
            if os.path.exists(fn1):    # @todo this catches files and directories
                return fn1
            # weird, typo is need to check the code tree
        url = url + '/example_data/' + fn
        print("url:",url)
        filename = url.split('/')[-1]
        if not os.path.exists(filename):    
            print(f"Downloading {filename}")
            wget.download(url,out=filename)
            print(f"\nRetrieved {filename}")
        else:
            print(f"{filename} already downloaded")
        return filename


    if sdfits != None:
        if dysh_data != None:
            fn = dysh_data + '/sdfits/' + sdfits
            if os.path.exists(fn):
                return fn
        print(f"could not handle sdfits={sdfits} yet")

    print("You have not given one of:   test=, example=, sdfits=")
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
