import inspect
import logging
import logging.config
import os
import sys
import time
from abc import ABC
from datetime import datetime
from functools import wraps
from io import StringIO
from pathlib import Path
from typing import Union

import _testcapi
from astropy.logger import AstropyLogger

try:
    DYSH_LOG_DIR = Path(os.getenv("DYSH_LOG_DIR", "."))
except ValueError as error:
    raise ValueError(f"DYSH_LOG_DIR must be set to a valid path! Got {DYSH_LOG_DIR}") from error

try:
    DYSH_HOME = Path(os.getenv("DYSH_HOME", "."))
except ValueError as error:
    raise ValueError(f"DYSH_HOME must be set to a valid path! Got {DYSH_HOME}") from error


class DatestampFileHandler(logging.FileHandler):
    def __init__(self, filename, header=None, mode="w", encoding=None, delay=0):
        if not filename:
            raise ValueError("Must provide a filename!")

        super().__init__(datetime.now().strftime(str(filename)), mode, encoding, delay)


class StringHandler(logging.StreamHandler):
    def __init__(self):
        string_stream = StringIO()
        print("initing string hamdler")
        super().__init__(stream=string_stream)


dysh_date_format = "%Y-%m-%dT%H:%M:%S%z"
dysh_history_formatter = logging.Formatter(
    fmt="%(asctime)s - %(levelname)s - %(modName)s.%(fName)s(%(message)s, %(extra)s)", datefmt=dysh_date_format
)


config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "extra_verbose": {
            "format": "%(asctime)s - %(pathname)s:%(lineno)s - %(message)s",
            "datefmt": dysh_date_format,
        },
        "verbose": {
            "format": "%(asctime)s - %(levelname)s - %(module)s - %(message)s",
            "datefmt": dysh_date_format,
        },
        "simple": {
            "format": "%(levelname)s %(message)s",
        },
        "super_simple": {
            "format": "%(message)s",
        },
        #        "data_history": {
        #            "format": "%(asctime)s - %(levelname)s - %(modName)s.%(fName)s(%(message)s, %(extra)s)",
        #           "datefmt": dysh_date_format,
        #        },
    },
    "handlers": {
        "stderr": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "simple",
        },
        "dysh_global_log_file": {
            "level": "INFO",
            "formatter": "verbose",
            "class": "logging.handlers.RotatingFileHandler",
            "filename": str(DYSH_HOME / "dysh.log"),
            "mode": "a",
        },
        "dysh_instance_log_file": {
            "level": "DEBUG",
            "formatter": "verbose",
            "class": "dysh.log.DatestampFileHandler",
            # Pull the log path from the environment. If this isn't set, an error
            # will be thrown. To log to the current directory, set this to .
            "filename": str(DYSH_LOG_DIR / "dysh.%Y_%m_%d_%H_%M_%S.log"),
            # Don't create the file until messages are written to it
            "delay": True,
        },
        #        "dysh_instance_string": {
        #            "level": "INFO",
        #           "formatter": "data_history",
        #            "class": "dysh.log.StringHandler",
        #        },
    },
    "loggers": {
        "dysh": {
            "handlers": [],
            "level": "WARNING",
        },
        #        "dysh_history": {
        #            "handlers": [],
        #            "level": "INFO",
        #        },
    },
}


def config_logging(verbosity: int, level: Union[int, None] = None, path: Union[Path, None] = None):
    global logger
    if verbosity:
        if verbosity == 0:
            level = logging.ERROR
        elif verbosity == 1:
            level = logging.WARNING
        elif verbosity == 2:
            level = logging.INFO
        elif verbosity == 3:
            level = logging.DEBUG
        else:
            raise ValueError(f"Invalid verbosity: {verbosity}")

    if level is None:
        raise ValueError("One of verbosity or level must be given!")

    config["loggers"]["dysh"]["handlers"] = ["stderr", "dysh_global_log_file"]
    # config["loggers"]["dysh_history"]["handlers"] = ["dysh_instance_string"]

    if path:
        config["handlers"]["dysh_instance_log_file"]["filename"] = path
        config["loggers"]["dysh"]["handlers"].append("dysh_instance_log_file")
        logger.info(f"Log file for this instance of dysh: {path}")

    # Init logging
    logging.config.dictConfig(config)
    logging.getLogger().setLevel(level)
    logger.setLevel(level)
    logger.debug(f"Logging has been set to verbosity {verbosity} / level {logging.getLevelName(level)}")


def log_function_call(func):
    """
    Decorator to log a function call

    Parameters
    ----------
    func : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except:  # remove the wrapper from the stack trace
            tp, exc, tb = sys.exc_info()
            _testcapi.set_exc_info(tp, exc, tb.tb_next)
            del tp, exc, tb
            raise
        # Log the function name and arguments
        sig = inspect.signature(func)
        logmsg = f"call: {func.__module__}"
        if hasattr(func, "__self__"):
            # func is  method of a class
            logmsg += f"{func.__self__.__class__.__name__}"
        logmsg += f".{func.__name__}{args}"
        if "kwargs" in sig.parameters:
            for k, v in kwargs.items():
                logmsg += f"{k}={v},"
        logger.info(logmsg)
        return result

    return wrapper


def format_dysh_log_record(record: logging.LogRecord) -> str:
    """
    Function to format a LogRecord into a string suitable
    for a FITS history card.  This function is intended for
    records of a dysh class method invocation and is used
    by :meth:`log_call_to_history`

    Parameters
    ----------
    record : logging.LogRecord
        The LogRecord to convert.   It will have some special
        purpose added attributes

    Returns
    -------
    str
        The history string

    """
    asctime = time.strftime(dysh_date_format, time.localtime(record.created))
    logmsg = f"{asctime} - {record.levelname} - {record.msg} : {record.modName}."
    if hasattr(record, "className"):
        logmsg += f"{record.className}."
    logmsg += f"{record.fName}("
    if len(record.args) > 0:
        logmsg += f"{record.args}"
    if len(record.kwargs) > 0:
        for k, v in record.kwargs.items():
            logmsg += f"{k}={v},"
    logmsg += ")"
    return logmsg


def log_call_to_history(func):
    """
    Decorator to log a class method call to the class's `_history` attribute. If the class
    has no such attribute, the function is still called but no logging takes place.

    Parameters
    ----------
    func : method
        The class method.

    Returns
    -------
    Any
        The result of the method call

    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self is None:  # not a class, but a function
            print(
                f"Function {func.__module__}.{func.__name__} is a function with no _history attribute. Use @log_function_call instead."
            )
            logger.warn(
                f"Function {func.__module__}.{func.__name__} is a function with no _history attribute. Use @log_function_call instead."
            )
            # could put this try around the whole thing but then it would catch exceptions from the wrapper itself
            try:
                result = func(*args, **kwargs)
            except:  # remove the wrapper from the stack trace
                tp, exc, tb = sys.exc_info()
                _testcapi.set_exc_info(tp, exc, tb.tb_next)
                del tp, exc, tb
                raise
        else:  # it's a class instance
            try:
                result = func(self, *args, **kwargs)
            except:  # remove the wrapper from the stack trace
                tp, exc, tb = sys.exc_info()
                _testcapi.set_exc_info(tp, exc, tb.tb_next)
                del tp, exc, tb
                raise
            classname = self.__class__.__name__
            print(f"got self {self} class {self.__class__} classname {classname}")
            if hasattr(self, "_history"):
                sig = inspect.signature(func)
                if "kwargs" in sig.parameters:
                    extra = {
                        "modName": func.__module__,
                        "className": classname,
                        "fName": func.__name__,
                        "kwargs": kwargs,
                    }
                else:
                    extra = {"modName": func.__module__, "className": classname, "fName": func.__name__, "extra": {}}
                with dhlogger.log_to_list() as log_list:
                    # dhlogger.handlers[0].setFormatter(dysh_history_formatter)
                    dhlogger.info("command", *args, extra=extra)
                    log_str = [format_dysh_log_record(i) for i in log_list]
                    self._history.extend(log_str)
            else:
                print(f"Class {func.__module__}.{classname} has no _history attribute. Use @log_function_call instead.")
                logger.warn(
                    f"Class {func.__module__}.{classname} has no _history attribute. Use @log_function_call instead."
                )
        return result

    return wrapper


class HistoricalBase(ABC):
    """Abstract base class to manage history and comments metadata."""

    def __init__(self):
        self._history = []
        self._comments = []
        print("Inst HistoricalBase")

    @property
    def history(self):
        """
        Get the history strings. These are typically converted to FITS HISTORY cards by the derived class.

        Returns
        -------
        list
            The list of string history

        """
        return self._history

    @property
    def comments(self):
        """
        Get the comment strings. These are typically converted to FITS COMMENT cards by the derived class.

        Returns
        -------
        list
            The list of string comments

        """
        return self._comments

    def add_comment(self, comment: str):
        """
        Add one or more comments to the class metadata.

        Parameters
        ----------
        comment : str or list of str
            The comment(s) to add

        Returns
        -------
        None.

        """
        if isinstance(comment, list):
            self._comments.extend(comment)
        else:
            self._comments.append(comment)

    def add_history(self, history: str):
        """
        Add one or more history entries to the class metadata

        Parameters
        ----------
        history : str or list of str
            The history card(s) to add

        Returns
        -------
        None.

        """
        # @todo should we prepend a datetime str here?
        if isinstance(history, list):
            self._history.extend(history)
        else:
            self._history.append(history)


logger = logging.getLogger("dysh")
dhlogger = AstropyLogger("dysh_history", level=logging.INFO)
