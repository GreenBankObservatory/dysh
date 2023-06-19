#!/usr/bin/env python
# NOTE: If any required directories are added, put them in MANIFEST.in or
# readthedocs build will break

from setuptools import setup, find_packages
import sys
import pathlib

def check_python(major,minor):
    try:
        assert sys.version_info >= (major,minor)
    except AssertionError:
        raise Exception("dysh requires you use Python %d.%d or above"%(major,minor))

# Ensure they are using Python 3.8 or above
check_python(3,8)

#excludelist= ["build","dist"]
excludelist= []
#print("Found packages ",find_packages(exclude=excludelist))

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="dysh",
    version = "0.1.0a3",
    author  = "Marc Pound",
    author_email = "mpound@umd.edu",
#    description = pdrtpy.DESCRIPTION,
#    keywords = pdrtpy.KEYWORDS,
     description="Dysh is a Python spectral line data reduction and analysis program for singledish data with specific emphasis on data from the Green Bank Telescope.",
    long_description=long_description, 
    long_description_content_type="text/markdown",
    packages = find_packages(exclude=excludelist),
    include_package_data = True,
    install_requires = [
        'astropy>=5.2.1',
        'numpy>=1.22.0',
        'scipy>=1.8.0',
        'matplotlib>=3.5.1',
        'pandas',
        'specutils',
    ],
    url = "https://pypi.org/project/dysh/",
    project_urls = {
        "Documentation": "https://dysh.readthedocs.io",
        "Source Code": "https://github.com/GreenBankObservatory/dysh",
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Intended Audience :: Science/Research",
    ],
    license = "GPLv3",
    zip_safe = False,
    python_requires = '>=3.8'
)
