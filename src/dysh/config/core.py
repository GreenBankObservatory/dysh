"""Configuration core functions"""

import contextlib
import importlib
import io
import os
import pkgutil
import warnings
from contextlib import nullcontext
from pathlib import Path
from textwrap import TextWrapper

from astropy import config as ac
from astropy.config import configuration as acc
from astropy.utils import silence
from astropy.utils.exceptions import AstropyDeprecationWarning

__all__ = ["create_config_file", "get_config_dir_path"]


def get_config_dir_path(rootname: str = "dysh") -> Path:
    """
    Determines the package configuration directory name and creates the
    directory if it doesn't exist.

    This directory is typically ``$HOME/.astropy/config``, but if the
    XDG_CONFIG_HOME environment variable is set and the
    ``$XDG_CONFIG_HOME/astropy`` directory exists, it will be that directory.
    If neither exists, the former will be created and symlinked to the latter.

    Parameters
    ----------
    rootname : str
        Name of the root configuration directory. For example, if ``rootname =
        'pkgname'``, the configuration directory would be ``<home>/.pkgname/``
        rather than ``<home>/.dysh`` (depending on platform).

    Returns
    -------
    configdir : Path
        The absolute path to the configuration directory.
    """

    return Path(ac.get_config_dir(rootname))


def recursive_subclasses(klass):
    """
    Yield all subclasses of the given class, per:
    https://adamj.eu/tech/2024/05/10/python-all-subclasses/
    """
    seen = set()
    for subclass in klass.__subclasses__():
        yield subclass
        for subsubclass in recursive_subclasses(subclass):
            if subsubclass not in seen:
                seen.add(subsubclass)
                yield subsubclass


def generate_config(pkgname="dysh", filename=None, verbose=False):
    """Generates a configuration file, from the list of `ConfigItem`
    objects for each subpackage.

    Parameters
    ----------
    pkgname : str or None
        The package for which to retrieve the configuration object.
    filename : str or file-like or None
        If None, the default configuration path is taken from `get_config`.

    """
    if verbose:
        verbosity = nullcontext
        filter_warnings = AstropyDeprecationWarning
    else:
        verbosity = silence
        filter_warnings = Warning

    package = importlib.import_module(pkgname)
    with verbosity(), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=filter_warnings)
        for mod in pkgutil.walk_packages(path=package.__path__, prefix=package.__name__ + "."):
            if mod.module_finder.path.endswith(("test", "tests")) or mod.name.endswith("setup_package"):
                # Skip test and setup_package modules
                continue
            if mod.name.split(".")[-1].startswith("_"):
                # Skip private modules
                continue

            with contextlib.suppress(ImportError):
                importlib.import_module(mod.name)

    wrapper = TextWrapper(initial_indent="## ", subsequent_indent="## ", width=78)

    if filename is None:
        filename = acc.get_config_filename(pkgname)

    with contextlib.ExitStack() as stack:
        if isinstance(filename, (str, os.PathLike)):
            fp = stack.enter_context(open(filename, "w"))
        else:
            # assume it's a file object, or io.StringIO
            fp = filename

        # Parse the subclasses, ordered by their module name
        # subclasses = ConfigNamespace.__subclasses__()
        subclasses = recursive_subclasses(ac.ConfigNamespace)
        processed = set()

        for conf in sorted(subclasses, key=lambda x: x.__module__):
            mod = conf.__module__

            # Skip modules for other packages, e.g. astropy modules that
            # would be imported when running the function for astroquery.
            if mod.split(".")[0] != pkgname:
                continue

            # Check that modules are not processed twice, which can happen
            # when they are imported in another module.
            if mod in processed:
                continue
            else:
                processed.add(mod)

            print_module = True
            for item in conf().values():
                if print_module:
                    # If this is the first item of the module, we print the
                    # module name, but not if this is the root package...
                    if item.module != pkgname:
                        modname = item.module.replace(f"{pkgname}.", "")
                        fp.write(f"[{modname}]\n\n")
                    print_module = False

                fp.write(wrapper.fill(item.description) + "\n")
                if isinstance(item.defaultvalue, (tuple, list)):
                    if len(item.defaultvalue) == 0:
                        fp.write(f"# {item.name} = ,\n\n")
                    elif len(item.defaultvalue) == 1:
                        fp.write(f"# {item.name} = {item.defaultvalue[0]},\n\n")
                    else:
                        fp.write(f"# {item.name} = {','.join(map(str, item.defaultvalue))}\n\n")
                else:
                    fp.write(f"# {item.name} = {item.defaultvalue}\n\n")


def create_config_file(pkg: str = "dysh", rootname: str = "dysh", overwrite: bool = False) -> bool:
    """
    Create the default configuration file for the specified package.
    If the file already exists, it is updated only if it has not been
    modified.  Otherwise the ``overwrite`` flag is needed to overwrite it.

    Parameters
    ----------
    pkg : str
        The package to be updated.
    rootname : str
        Name of the root configuration directory.
    overwrite : bool
        Force updating the file if it already exists.

    Returns
    -------
    updated : bool
        If the profile was updated, `True`, otherwise `False`.

    """
    # local import to prevent using the logger before it is configured
    from astropy.logger import log

    cfgfn = Path(acc.get_config_filename(pkg, rootname=rootname))

    # generate the default config template
    template_content = io.StringIO()
    generate_config(pkg, template_content)
    template_content.seek(0)
    template_content = template_content.read()

    doupdate = True

    # if the file already exists, check that it has not been modified
    if cfgfn is not None and cfgfn.is_file():
        with open(cfgfn, encoding="latin-1") as fd:
            content = fd.read()

        doupdate = acc.is_unedited_config_file(content, template_content)

    if doupdate or overwrite:
        with open(cfgfn, "w", encoding="latin-1") as fw:
            fw.write(template_content)
        log.info(f"The configuration file has been successfully written to {cfgfn}")
        return True
    elif not doupdate:
        log.warning(
            "The configuration file already exists and seems to "
            "have been customized, so it has not been updated. "
            "Use overwrite=True if you really want to update it."
        )

    return False
