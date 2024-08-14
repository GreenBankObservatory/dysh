import inspect
import logging
import logging.config
import os
from datetime import datetime
from functools import wraps
from io import StringIO
from pathlib import Path
from typing import Union

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


##class DyshHistoryLogger():
#  def __init
config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "extra_verbose": {
            "format": "%(asctime)s - %(pathname)s:%(lineno)s - %(message)s",
            "datefmt": "%Y-%m-%dT%H:%M:%S%z",
        },
        "verbose": {
            "format": "%(asctime)s - %(levelname)s - %(module)s - %(message)s",
            "datefmt": "%Y-%m-%dT%H:%M:%S%z",
        },
        "simple": {
            "format": "%(levelname)s %(message)s",
        },
        "super_simple": {
            "format": "%(message)s",
        },
        "data_history": {
            "format": "%(asctime)s - %(levelname)s - %(modName)s.%(fName)s(%(message)s, %(extra)s)",
            "datefmt": "%Y-%m-%dT%H:%M:%S%z",
        },
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
        "dysh_instance_string": {
            "level": "INFO",
            "formatter": "data_history",
            "class": "dysh.log.StringHandler",
        },
    },
    "loggers": {
        "dysh": {
            "handlers": [],
            "level": "WARNING",
        },
        "dysh_history": {
            "handlers": [],
            "level": "INFO",
        },
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
    config["loggers"]["dysh_history"]["handlers"] = ["dysh_instance_string"]

    if path:
        config["handlers"]["dysh_instance_log_file"]["filename"] = path
        config["loggers"]["dysh"]["handlers"].append("dysh_instance_log_file")
        logger.info(f"Log file for this instance of dysh: {path}")

    # Init logging
    logging.config.dictConfig(config)
    logging.getLogger().setLevel(level)
    logger.setLevel(level)
    logger.debug(f"Logging has been set to verbosity {verbosity} / level {logging.getLevelName(level)}")
    dhlogger.setLevel(logging.INFO)


def log_call(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        # Log the function name and arguments
        # dhlogger.info(f"calling {args}", extra=kwargs)
        sig = inspect.signature(func)
        print(sig)
        if "kwargs" in sig.parameters:
            extra = {"modName": func.__module__, "fName": func.__name__, "extra": kwargs}
        else:
            extra = {"modName": func.__module__, "fName": func.__name__, "extra": ""}
        dhlogger.info(*args, extra=extra)
        # Call the original function
        result = func(*args, **kwargs)
        # Return the result
        return result

    return wrap


logger = logging.getLogger("dysh")
dhlogger = logging.getLogger("dysh_history")
