# Dysh

Dysh is a Python spectral line data reduction and analysis program for singledish data with specific emphasis on data from the Green Bank Telescope.  It is currently under development in collaboration between the Green Bank Observatory and the Laboratory for Millimeter-Wave Astronomy (LMA) at University of Maryland (UMD).  It is intended to be a full replacement for the GBO's current reduction package GBTIDL.

## Installation

Dysh requires Python 3.8+ and recent versions of [astropy]( https://astropy.org), [numpy](https://numpy.org), [scipy](https://scipy.org), [pandas](https://pandas.pydata.org) and [matplotlib](https://matplotlib.org). Once you are in a suitable environment, simply:

```bash
    $ pip install dysh
```


## Getting Help & Giving Feedback

more here

## Reporting Issues

If you find a bug or something you think is in error, please report it on
the [github issue tracker](https://github.com/GreenBankObservatory/dysh/issues).
(You must have a [Github account](https://github.com) to submit an issue)

---

## Development

Here are the steps if you want to develop code for dysh.  We use [hatch](https://hatch.pypa.io/) to manage the build environment.

1.  Clone the repo and install hatch.

```bash
    $ git clone git@github.com:GreenBankObservatory/dysh.git
    $ cd dysh
    $ pip install hatch  # or pipx if you prefer
```

2.  Create and activate a virtual environment with hatch and install the packages required for development.
The virtual environment will be created the first time; subsequent invoking ``hatch shell`` will simply load the created environment.cdi

```bash
    $ hatch shell
    $ pip install -r requirements_dev.txt
```

3.  Build and install the package

```bash
    $ hatch build
    $ pip install -e .
```


#### pre-commit

This repository provides pre-commit hooks that enforce some code formatting/quality checks. To use them:

1. Install [pre-commit](https://pre-commit.com/)
2. Install the hooks via <font size="5">`pre-commit install`</font>
3. Future commits will run these hooks prior to allowing the commit

## Testing

```bash
    $ pytest
```

