[![Documentation Status](https://readthedocs.org/projects/dysh/badge/?version=latest)](https://dysh.readthedocs.io/en/latest/?badge=latest)
[![pre-commit.ci Status](https://results.pre-commit.ci/badge/github/GreenBankObservatory/dysh/main.svg)](https://results.pre-commit.ci/latest/github/GreenBankObservatory/dysh/main)
[![CI Workflow Build Status](https://github.com/GreenBankObservatory/dysh/actions/workflows/ci.yml/badge.svg)](https://github.com/GreenBankObservatory/dysh/actions/workflows/ci.yml)

# dysh

*dysh* is a Python spectral line data reduction and analysis program for singledish data with specific emphasis on data from the Green Bank Telescope. It is currently under development in collaboration between the [Green Bank Observatory](https:/greenbankobservatory.org) and the Laboratory for Millimeter-Wave Astronomy (LMA) at [University of Maryland (UMD)](https://www.astro.umd.edu). It is intended to be an alternative to the GBO's current reduction package [GBTIDL](https://www.gb.nrao.edu/GBT/DA/gbtidl/users_guide/).

## Installation

Note: if you are on the GBO network, dysh will already be installed; you do not need to do anything further! Other uses cases are outlined below.

### Global Installation via pipx

Example use case: you want to quickly install and use dysh on a non-GBO computer

If you want to avoid managing a virtual environment, the easiest way is via [pipx](https://github.com/pypa/pipx#install-pipx). After installing pipx, dysh can be installed via:

```sh
# Install dysh
$ pipx install dysh[nb]
# Launch dysh notebook
$ dysh-lab
```

### Local Installation via pip

Example use case: you want to use dysh with a specific set of dependencies, or in conjunction with an existing project

#### Stable Version
dysh is most easily installed with *pip*, which will take care of any dependencies. The packaged code is hosted at the [Python Packaging Index](https://pypi.org/project/dysh).

```sh
$ pip install dysh
```

#### Beta Version

Beta versions will also be published to PyPI, and can be installed via:

```sh
$ pip install dysh --pre
```

#### Development Version

Development versions can be installed from GitHub branches via:

```sh
$ pip install git+https://github.com/GreenBankObservatory/dysh
```
For more options, see the [pip VCS Support documentation](https://pip.pypa.io/en/stable/topics/vcs-support/).

## Reporting Issues

If you find a bug or something you think is in error, please report it on
the [GitHub issue tracker](https://github.com/GreenBankObservatory/dysh/issues).
(You must have a [GitHub account](https://github.com) to submit an issue)

---

## Development

See the [For Developers](https://dysh.readthedocs.io/en/latest/for_developers/index.html) documentation for more detailed instructions on setting up a development environment.

### Clone the Repo

```sh
$ git clone git@github.com:GreenBankObservatory/dysh.git
```

### Set Up Dysh Environment

#### Via Hatch

The recommended development workflow is to use [hatch](https://hatch.pypa.io/latest/tutorials/environment/basic-usage/). After installing hatch, this will look something like:

```sh
$ hatch shell
```

#### Without Hatch

If you do not want to use Hatch, it is possible to develop using a "classic" workflow. From the root of the dysh repo:

```sh
$ # Create your virtual environment
$ python -m venv /path/to/venv
$ # Activate your virtual environment
$ source /path/to/venv/bin/activate
$ # Install dysh and its development dependencies
$ pip install -e .[all]
```

### Testing
We use pytest for unit and integration testing. From the top-level dysh directory, run:

```sh
$ pytest
```
