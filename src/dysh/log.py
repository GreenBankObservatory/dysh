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
from typing import Callable, NewType, Union  # , Self # not available until 3.11

from astropy import log
from astropy.io.fits.header import _HeaderCommentaryCards
from astropy.logger import AstropyLogger

from . import version
from .ascii import ensure_ascii

# Set Astropy logging level to warning.
log.setLevel("WARNING")

LOGGING_INITIALIZED = False
logger = logging.getLogger("dysh")
dhlogger = AstropyLogger("dysh_history", level=logging.INFO)

_DYSH_LOG_DIR = os.getenv("DYSH_LOG_DIR", ".")
try:
    DYSH_LOG_DIR = Path(_DYSH_LOG_DIR)
except ValueError as error:
    raise ValueError(f"DYSH_LOG_DIR must be set to a valid path! Got {_DYSH_LOG_DIR}") from error

_DYSH_HOME = os.getenv("DYSH_HOME", ".")
try:
    DYSH_HOME = Path(_DYSH_HOME)
except ValueError as error:
    raise ValueError(f"DYSH_HOME must be set to a valid path! Got {_DYSH_HOME}") from error


class DatestampFileHandler(logging.FileHandler):
    def __init__(self, filename, mode="w", encoding=None, delay=False):
        if not filename:
            raise ValueError("Must provide a filename!")

        super().__init__(filename=datetime.now().strftime(str(filename)), mode=mode, encoding=encoding, delay=delay)


class StringHandler(logging.StreamHandler):
    def __init__(self):
        string_stream = StringIO()
        super().__init__(stream=string_stream)


dysh_date_format = "%Y-%m-%dT%H:%M:%S%z"
dysh_version = version()
dysh_history_formatter = logging.Formatter(
    fmt="%(asctime)s - %(levelname)s - %(modName)s.%(fName)s(%(message)s, %(extra)s)", datefmt=dysh_date_format
)


def dysh_date():
    return datetime.now().strftime(dysh_date_format)


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


def init_logging(verbosity: int, level: Union[int, None] = None, path: Union[Path, None] = None, quiet=False):
    global LOGGING_INITIALIZED
    if LOGGING_INITIALIZED is True:
        logger.warning(
            "dysh logging has already been initialized! Continuing, but this behavior is not well defined, "
            "and you will likely end up with duplicate log handlers"
        )

    LOGGING_INITIALIZED = True
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

    if quiet:
        config["handlers"]["stderr"]["level"] = "WARNING"

    if path:
        config["handlers"]["dysh_instance_log_file"]["filename"] = path
        config["loggers"]["dysh"]["handlers"].append("dysh_instance_log_file")
        logger.info(f"Log file for this instance of dysh: {path}")

    # Init logging
    logging.config.dictConfig(config)
    logging.getLogger().setLevel(level)
    logger.setLevel(level)
    logger.debug(f"Logging has been set to verbosity {verbosity} / level {logging.getLevelName(level)}")


def log_function_call(log_level: str = "info"):
    """
    Decorator to log a function call

    Parameters
    ----------
    log_level : str
        The logging level to use for logging. One of
        ['CRITICAL', 'FATAL', 'ERROR', 'WARN', 'WARNING',
        'INFO', 'DEBUG', 'NOTSET'].
        Case-insensitive. Default: 'info'

    Returns
    -------
    inner_decorator : Any
        Whatever the function returns.

    """
    # if not callable(arg):  doesn't work
    #    log_level = "info"

    def inner_decorator(func: Callable):
        # the inner decorator is to process the log_level argument of
        # the outer decorator

        try:
            ilog_level = logging._nameToLevel[log_level.upper()]
        except KeyError:
            raise Exception(f"Log level {log_level} unrecognized. Must be one of {list(logging._nameToLevel.keys())}.")  # noqa: B904

        @wraps(func)
        def func_wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
            except:  # remove the wrapper from the stack trace  # noqa: E722
                tp, exc, tb = sys.exc_info()
                raise tp(exc).with_traceback(tb.tb_next)  # noqa: B904
            # Log the function name and arguments
            sig = inspect.signature(func)
            logmsg = f"DYSH v{dysh_version} : {func.__module__}"
            if hasattr(func, "__self__"):
                # func is  method of a class
                logmsg += f"{func.__self__.__class__.__name__}"
            logmsg += f".{func.__name__}{args}"
            if "kwargs" in sig.parameters:
                for k, v in kwargs.items():
                    logmsg += f"{k}={v},"
            logmsg = ensure_ascii(logmsg)
            logger.log(level=ilog_level, msg=logmsg)
            return result

        return func_wrapper

    return inner_decorator


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
    # NB: possibly get rid of levelname.  Also possibly move dysh_version to before asctime.
    # logmsg = f"{asctime} - {record.levelname} - {record.msg} : {record.modName}."
    logmsg = f"{asctime} - {record.msg} : {record.modName}."
    if hasattr(record, "className"):
        logmsg += f"{record.className}."
    logmsg += f"{record.fName}("
    if len(record.args) > 0:
        for a in record.args:
            logmsg += f"{a},"
    if hasattr(record, "kwargs"):
        if len(record.kwargs) > 0:
            for k, v in record.kwargs.items():
                logmsg += f"{k}={v},"
    logmsg += ")"
    return ensure_ascii(logmsg)


def log_call_to_result(func: Callable):
    """
    Decorator to log a class method call to the method's resultant class's `_history` attribute.
    If the resultant class has no such attribute, the function is still called but no logging takes place.

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
        if self is None:
            try:
                result = func(*args, **kwargs)
            except:  # remove the wrapper from the stack trace  # noqa: E722
                # @todo this no longer works under python >3.9
                tp, exc, tb = sys.exc_info()
                raise tp(exc).with_traceback(tb.tb_next)  # noqa: B904
        else:
            try:
                result = func(self, *args, **kwargs)
            except:  # remove the wrapper from the stack trace  # noqa: E722
                tp, exc, tb = sys.exc_info()
                raise tp(exc).with_traceback(tb.tb_next)  # noqa: B904
        resultname = result.__class__.__name__
        if hasattr(result, "_history"):
            sig = inspect.signature(func)
            if self is not None:
                extra = {
                    "modName": func.__module__,
                    "className": self.__class__.__name__,
                    "fName": func.__name__,
                }
            else:
                extra = {
                    "modName": func.__module__,
                    "fName": func.__name__,
                }
            if "kwargs" in sig.parameters:
                extra["kwargs"] = kwargs
            with dhlogger.log_to_list() as log_list:
                dhlogger.info(f"DYSH v{dysh_version}", *args, extra=extra)
                log_str = [format_dysh_log_record(i) for i in log_list]
                if hasattr(result, "_history"):
                    result._history.extend(log_str)
        else:
            logger.warn(f"Class {resultname} has no _history attribute. Use @log_function_call instead.")
        return result

    return wrapper


def log_call_to_history(func: Callable):
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
            logger.warn(
                f"Function {func.__module__}.{func.__name__} is a function with no _history attribute. Use"
                " @log_function_call instead."
            )
            # could put this try around the whole thing but then it would catch exceptions from the wrapper itself
            try:
                result = func(*args, **kwargs)
            except:  # remove the wrapper from the stack trace  # noqa: E722
                tp, exc, tb = sys.exc_info()
                raise tp(exc).with_traceback(tb.tb_next)  # noqa: B904
        else:  # it's a class instance
            try:
                result = func(self, *args, **kwargs)
            except:  # remove the wrapper from the stack trace  # noqa: E722
                tp, exc, tb = sys.exc_info()
                raise tp(exc).with_traceback(tb.tb_next)  # noqa: B904
            classname = self.__class__.__name__
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
                    dhlogger.info(f"DYSH v{dysh_version}", *args, extra=extra)
                    log_str = [format_dysh_log_record(i) for i in log_list]
                    self._history.extend(log_str)
                    # Handle the case where we want the result, e.g. a ScanBlock
                    # to record the parent class call that created it.
                    if hasattr(result, "_history"):
                        result._history.extend(log_str)
            else:
                logger.warn(
                    f"Class {func.__module__}.{classname} has no _history attribute. Use @log_function_call instead."
                )
        return result

    return wrapper


# type History = list[str]  # not available till python 3.12
StrList = NewType("Strlist", list[str])


class HistoricalBase(ABC):  # noqa: B024
    """Abstract base class to manage history and comments metadata."""

    def __init__(self):
        self._history = StrList([])
        self._comments = StrList([])

    def _remove_duplicates(self):
        # Use dict.fromkeys here to preserve order. set() does not.
        self._history = list(dict.fromkeys(self._history))
        self._comments = list(dict.fromkeys(self._comments))

    @property
    def history(self) -> StrList:
        """
        Get the history strings. These are typically converted to FITS HISTORY cards by the derived class.
        Duplicate entries are removed.  History strings are cleaned of non-ASCII and
        non-printable characters that are invalid in FITS Cards

        Returns
        -------
        list
            The list of string history

        """
        # remove duplicates, due to inherited classes
        self._remove_duplicates()
        # time sort
        return sorted(ensure_ascii(self._history))

    @property
    def comments(self) -> StrList:
        """
        Get the comment strings. These are typically converted to FITS COMMENT cards by the derived class.
        Duplicate comments are removed. Comment  strings are cleaned of non-ASCII and
        non-printable characters that are invalid in FITS Cards

        Returns
        -------
        list
            The list of string comments

        """
        # remove duplicates, due to inherited classes
        self._remove_duplicates()
        return ensure_ascii(self._comments)

    def add_comment(self, comment: Union[str, StrList], add_time: bool = False) -> None:
        """
        Add one or more comments to the class metadata.

        Parameters
        ----------
        comment : str or list of str or `~astropy.io.fits.header._HeaderCommentaryCards`
            The comment(s) to add
        add_time: bool
            If True, prepend the date and time the history was added

        Returns
        -------
        None.

        """
        if add_time:
            _time = datetime.now().strftime(dysh_date_format)
        if isinstance(comment, list) or isinstance(comment, _HeaderCommentaryCards):
            if add_time:  # tag each comment in list
                for c in comment:
                    h = f"{_time} - {c}"
                    self.add_comment(h, add_time=False)
            else:
                self._comments.extend(comment)
        else:
            if add_time:
                h = f"{_time} - {comment}"
                self._comments.append(h)
            else:
                self._comments.append(comment)

    def add_history(self, history: Union[str, StrList], add_time: bool = False) -> None:
        """
        Add one or more history entries to the class metadata

        Parameters
        ----------
        history : str or list of str or `~astropy.io.fits.header._HeaderCommentaryCards`
            The history card(s) to add
        add_time: bool
            If True, prepend the date and time the history was added.

        Returns
        -------
        None.

        """
        history = ensure_ascii(history, check=True)
        if add_time:
            _time = datetime.now().strftime(dysh_date_format)
        if isinstance(history, list) or isinstance(history, _HeaderCommentaryCards):
            for c in history:
                if add_time:  # tag each comment in list
                    h = f"{_time} - {c}"
                    self.add_history(h, add_time=False)
                else:
                    self._history.extend(history)
        else:
            if add_time:
                h = f"{_time} - {history}"
                self._history.append(h)
            else:
                self._history.append(history)

    # def merge_commentary(self, other: Self) -> None:  # Self not available until python 3.11
    def merge_commentary(self, other: object) -> None:
        """
        Merge the history and comments from another HistoricalBase instance.
        The history and comments are added to this instance and duplicates
        are removed.

        Parameters
        ----------
        other : ~log.HistoricalBase
            An class instance that has history and comment attributes.

        Returns
        -------
        None

        """
        if hasattr(other, "_history"):
            self.add_history(other._history)
        if hasattr(other, "_comments"):
            self.add_comment(other._comments)
        self._remove_duplicates()
