# Dysh

Dysh will be a replacement for GBTIDL. It is currently under development in collaboration with The Laboratory for Millimeter-Wave Astronomy (LMA) at University of Maryland (UMD).

## Installation

Dysh requires Python 3.8+. Once you are in a suitable environment, simply:

```bash
$ pip install -e .
```


## Development

If you are working on Dysh itself, here's how to get your environment set up

First, you'll need `hatch` installed. The "proper" way to do that is via pipx:

```bash
$ pipx install hatch
```

But you can also manage your virtual environment entirely yourself:

```bash
$ python3.8 -m venv /path/to/venvs/dysh-3.8
$ source /path/to/venvs/dysh-3.8/bin/activate
$ pip install hatch
```

Once you have `hatch`, simply:

```
$ pip install -e .
```

To validate the install, you can run Dysh CLI via:


```bash
$ dysh --help
```

## Testing

Simply:

```bash
$ pytest
```

### pre-commit

This repository provides pre-commit hooks that enforce some code formatting/quality checks. To use them:

1. Install [pre-commit](https://pre-commit.com/)
2. Install the hooks via `pre-commit install`
3. Future commits will run these hooks prior to allowing the commit
