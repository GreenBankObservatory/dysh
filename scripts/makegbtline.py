#!/usr/bin/env python3
# make something like the GBTIDL RR_w_KuP file
import argparse
import sys

import astropy.units as u
import pandas as pd
from astropy.table import vstack
from astroquery.splatalogue import Splatalogue

from dysh.util import get_project_data, replace_col_astype

parser = argparse.ArgumentParser(
    prog=sys.argv[0],
    description="Create spectral line or recombination line CSV file for use with dysh's spectral line search feature.",
)
parser.add_argument(
    "--line", "-l", action="store", help="Output file name for spectral line data", default=None, required=False
)
parser.add_argument(
    "--recomb", "-r", action="store", help="Output file name for recombination line data", default=None, required=False
)
parser.add_argument(
    "--overwrite", "-w", action="store_true", help="overwrite a previous output file (astropy Table)", required=False
)
parser.add_argument("--verbose", "-v", action="store_true", help="verbose output", default=False, required=False)
parser.add_argument(
    "--redshift",
    "-z",
    action="store_true",
    help="add bright high frequency lines that could be redshifted into the GBT bands",
    default=False,
    required=False,
)

# this is the original line from grom GBTIDL. We use this to get the names of the molecules to search for.
gbtidl_file = get_project_data() / "GBTIDL_RRF_w_Kup.csv.gz"
# add $ to regex to eliminate, e.g., "Carbon Monoxide Ion" as a match
redshift_lines = [
    "Carbon Monoxide$",
    "Atomic Hydrogen",
    "Atomic Carbon",
    "Atomic Oxygen",
    "Atomic Helium",
    "Atomic Deuterium",
    "Carbon Monosulfide$",
    r"HCO\+",
    " HCN ",
]
redshift_recomb_lines = ["Recombination"]

args = parser.parse_args()
lowfreq = 289 * u.MHz
maxfreq = 120 * u.GHz

if args.line is None and args.recomb is None:
    print("No output file(s) specified.  Specify at least one of -l or -r.")
    exit(0)

if args.line is not None:
    p = pd.read_csv(
        gbtidl_file,
        delimiter=":",
        names=["rf", "rfe", "species", "name", "qn", "intensity", "eu", "intref", "freq"],
        comment=";",
    )

    x = list(set(p["name"]))
    # sx = list(set([s.split()[0] for s in x]))
    sx = list(set(x))
    skip = []
    print(f"Total species: {len(sx)}")
    t = Splatalogue.query_lines(
        min_frequency=lowfreq,
        max_frequency=maxfreq,
        chemical_name=f"{sx[0]}",
        intensity_lower_limit=-9,
        intensity_type="CDMS/JPL (log)",
    )
    for i in sx[1:]:
        try:
            if "Atomic" in i:
                ill = None
            else:
                ill = -9
            t2 = Splatalogue.query_lines(
                min_frequency=lowfreq,
                max_frequency=maxfreq,
                chemical_name=f"{i}",
                intensity_lower_limit=ill,
                intensity_type="CDMS/JPL (log)",
            )
        except Exception as e:
            if args.verbose:
                print(f"skipping {i} because {e}")
            skip.append(i)
            continue
        t = vstack([t, t2])

    if args.redshift:
        if args.verbose:
            print(f"Adding high-z transition for {redshift_lines}")
        for k in redshift_lines:
            t2 = Splatalogue.query_lines(
                maxfreq, 60 * u.micron, chemical_name=k, intensity_lower_limit=-9, intensity_type="CDMS/JPL (log)"
            )
            t = vstack([t, t2])

    # now replace intintensity string column with a float.
    replace_col_astype(t, "intintensity", float, -1e20)
    t.write(args.line, format="ascii.ecsv", overwrite=args.overwrite)
    print(
        f"Total species written {len(set(sx) - set(skip))}.\nNumber species skipped: {len(skip)}.\nLength of final line table {args.line} = {len(t)}"
    )
if args.recomb is not None:
    if args.verbose:
        print("Doing recomb lines")
    if args.redshift:
        maxfreq = 60 * u.micron
        if args.verbose:
            print("Including high-z recomb lines")
    t = Splatalogue.query_lines(lowfreq, maxfreq, chemical_name="Recombination")
    replace_col_astype(t, "intintensity", float, -1e20)
    t.write(args.recomb, format="ascii.ecsv", overwrite=args.overwrite)
    print(f"Length of final recomb table {args.recomb} = {len(t)}")
