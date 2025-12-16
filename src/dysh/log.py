import inspect
import logging
import logging.config
import os
import sys
import time
from abc import ABC
from collections.abc import Callable  # , Self # not available until 3.11
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import NewType

from astropy import log as astropy_logger
from astropy.io.fits.header import _HeaderCommentaryCards
from astropy.logger import AstropyLogger
from platformdirs import user_log_dir

from . import version
from .ascii import ensure_ascii

# Set Astropy logging level to warning.
astropy_logger.setLevel("WARNING")

DYSH = "dysh"
USER_LOG_DIR = Path(user_log_dir(DYSH))

# avoid "No handler" warnings in certain contexts
logging.getLogger(DYSH).addHandler(logging.NullHandler())

logger = logging.getLogger(DYSH)
logger._configured = False
dhlogger = AstropyLogger("dysh_history", level=logging.INFO)
dhlogger._configured = False

global_logger = logging.getLogger(f"{DYSH}_global")
instance_logger = logging.getLogger(f"{DYSH}_instance")
# Environment variable to disable logging decorators for performance benchmarking
DISABLE_HISTORY_LOGGING = os.getenv("DYSH_DISABLE_HISTORY_LOGGING", "0") == "1"


dysh_date_format = "%Y-%m-%dT%H:%M:%S%z"
FORMATTERS = {
    "verbose": logging.Formatter(
        "%(asctime)s %(levelname).1s %(process)d %(name)s:%(lineno)d %(message)s",
        "%Y-%m-%dT%H:%M:%S%z",
    ),
    "simple": logging.Formatter("%(asctime)s.%(msecs).3d %(levelname).1s %(message)s", "%H:%M:%S"),
    "history": logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(modName)s.%(fName)s(%(message)s, %(extra)s)", datefmt=dysh_date_format
    ),
}
dysh_version = version()


def dysh_date():
    """A date formatted for dysh log messages: '%Y-%m-%dT%H:%M:%S%z'"""
    return datetime.now().strftime(dysh_date_format)


def init_console_log(level):
    try:
        level_val = getattr(logging, str(level).upper())
    except AttributeError:
        level_val = int(level)

    logger.setLevel(level_val)

    console_handler = logging.StreamHandler(stream=sys.stderr)
    console_handler.setLevel(level_val)
    console_handler.setFormatter(FORMATTERS["simple"])
    logger.addHandler(console_handler)


def init_global_log(
    global_log_path: Path | None = None, rotate_bytes: int = 10_000_000, backups: int = 5, level: str | int = "INFO"
) -> Path:
    log_path = Path(global_log_path) if global_log_path else USER_LOG_DIR / f"{DYSH}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.handlers.RotatingFileHandler(
        log_path, maxBytes=rotate_bytes, backupCount=backups, encoding="utf-8"
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(FORMATTERS["verbose"])
    # Make sure we can actually write it!
    log_path.touch()
    logger.addHandler(file_handler)
    global_logger.addHandler(file_handler)
    global_logger.setLevel(logging.DEBUG)

    return log_path


def init_instance_log(
    instance_log_dir: Path = Path("."), instance_log_file: Path | None = None, level: str | int = "INFO"
) -> Path:
    instance_log_dir.mkdir(parents=True, exist_ok=True)
    if instance_log_file is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = instance_log_dir / f"{DYSH}-{os.getpid()}-{ts}.log"
    else:
        log_path = instance_log_dir / instance_log_file
    # Make sure we can actually write it!
    log_path.touch()
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(FORMATTERS["verbose"])
    logger.addHandler(file_handler)
    instance_logger.addHandler(file_handler)
    instance_logger.setLevel(logging.DEBUG)
    return log_path


def init_logging(verbosity: int | None = None, level: int | None = None, path: Path | None = None, quiet=False):
    # Clear existing handlers to avoid duplicates when re-initializing
    if logger._configured:
        logger.handlers.clear()
    if dhlogger._configured:
        dhlogger.handlers.clear()

    # If verbosity not provided, check environment variable
    if verbosity is None:
        verbosity = int(os.environ.get("DYSH_VERBOSITY", "2"))

    if quiet:
        level = logging.WARNING
    elif verbosity == 0:
        level = logging.ERROR
    elif verbosity == 1:
        level = logging.WARNING
    elif verbosity == 2:
        level = logging.INFO
    elif verbosity == 3:
        level = logging.DEBUG
    else:
        raise ValueError(f"Invalid verbosity: {verbosity}")

    init_console_log(level)
    logger.debug(f"Logging has been set to verbosity {verbosity} / level {logging.getLevelName(level)}")
    try:
        global_log_path = init_global_log()
    except Exception as e:
        logger.warning(f"Failed to initialize dysh_global_log_file: {e}")
    else:
        logger.debug(f"Global file for {DYSH}: {global_log_path}")

    if path:
        try:
            instance_log_path = init_instance_log(instance_log_dir=path.parent, instance_log_file=path.name)
        except Exception as e:
            raise ValueError(f"Failed to initialize requested instance log file in {path}") from e
        logger.info(f"Log file for this instance of {DYSH}: {instance_log_path.absolute()}")
    else:
        logger.debug("Instance log file not requested; disabling")

    logger._configured = True
    dhlogger._configured = True


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
        # If logging is disabled, return the function unchanged
        if DISABLE_HISTORY_LOGGING:
            return func

        # the inner decorator is to process the log_level argument of
        # the outer decorator

        try:
            ilog_level = logging._nameToLevel[log_level.upper()]
        except KeyError:
            raise Exception(f"Log level {log_level} unrecognized. Must be one of {list(logging._nameToLevel.keys())}.")  # noqa: B904

        # Cache signature at decoration time for performance
        sig = inspect.signature(func)

        @wraps(func)
        def func_wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
            except:  # noqa: E722 - intentionally catch all to manipulate traceback
                _tp, exc, tb = sys.exc_info()
                # Re-raise with modified traceback to hide this wrapper from stack trace
                raise exc.with_traceback(tb.tb_next) from None
            # Log the function name and arguments
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

    # If logging is disabled, return the function unchanged
    if DISABLE_HISTORY_LOGGING:
        return func

    # Cache signature at decoration time for performance
    sig = inspect.signature(func)

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self is None:
            try:
                result = func(*args, **kwargs)
            except:  # noqa: E722 - intentionally catch all to manipulate traceback
                _, exc, tb = sys.exc_info()
                # Re-raise with modified traceback to hide this wrapper from stack trace
                raise exc.with_traceback(tb.tb_next) from None
        else:
            try:
                result = func(self, *args, **kwargs)
            except:  # noqa: E722 - intentionally catch all to manipulate traceback
                _tp, exc, tb = sys.exc_info()
                # Re-raise with modified traceback to hide this wrapper from stack trace
                raise exc.with_traceback(tb.tb_next) from None
        resultname = result.__class__.__name__
        if hasattr(result, "_history"):
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
            logger.warning(f"Class {resultname} has no _history attribute. Use @log_function_call instead.")
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

    # If logging is disabled, return the function unchanged
    if DISABLE_HISTORY_LOGGING:
        return func

    # Cache signature at decoration time for performance
    sig = inspect.signature(func)

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self is None:  # not a class, but a function
            logger.warning(
                f"Function {func.__module__}.{func.__name__} is a function with no _history attribute. Use"
                " @log_function_call instead."
            )
            # could put this try around the whole thing but then it would catch exceptions from the wrapper itself
            try:
                result = func(*args, **kwargs)
            except:  # noqa: E722 - intentionally catch all to manipulate traceback
                _, exc, tb = sys.exc_info()
                # Re-raise with modified traceback to hide this wrapper from stack trace
                raise exc.with_traceback(tb.tb_next) from None
        else:  # it's a class instance
            try:
                result = func(self, *args, **kwargs)
            except:  # noqa: E722 - intentionally catch all to manipulate traceback
                _, exc, tb = sys.exc_info()
                # Re-raise with modified traceback to hide this wrapper from stack trace
                raise exc.with_traceback(tb.tb_next) from None
            classname = self.__class__.__name__
            if hasattr(self, "_history"):
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
                    dhlogger.info(f"DYSH v{dysh_version}", *args, extra=extra)
                    log_str = [format_dysh_log_record(i) for i in log_list]
                    self._history.extend(log_str)
                    # Handle the case where we want the result, e.g. a ScanBlock
                    # to record the parent class call that created it.
                    if hasattr(result, "_history"):
                        result._history.extend(log_str)
            else:
                logger.warning(
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

    def add_comment(self, comment: str | StrList, add_time: bool = False) -> None:
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
        elif add_time:
            h = f"{_time} - {comment}"
            self._comments.append(h)
        else:
            self._comments.append(comment)

    def add_history(self, history: str | StrList, add_time: bool = False) -> None:
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
        elif add_time:
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
