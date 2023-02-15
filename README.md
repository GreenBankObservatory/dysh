# astrodysh

astrodysh!

## Installation

astrodysh requires Python 3.11+. Once you are in a suitable environment, simply:

```bash
# Currently available only via GBO's private PyPI repository
$ pip install pyspeckit --extra-index-url http://pypi.gb.nrao.edu/simple
```

## Development

If you are working on astrodysh itself, here's how to get your environment set up

```bash
$ python3.11 -m venv /path/to/venvs/astrodysh-3.11
$ source /path/to/venvs/astrodysh-3.11/bin/activate
$ pip install -U pip setuptools wheel pdm
# pyspeckit must be manually installed via pip to avoid errors :(
$ pip install pyspeckit
$ pdm install
```

To validate your virtual environment, you can run the tests:

```bash
$ pytest
```

And run the CLI directly, via:

```bash
$ astrodysh --help
```
