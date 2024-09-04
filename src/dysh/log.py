import logging
import logging.config
import os
from datetime import datetime
from pathlib import Path
from typing import Union

LOGGING_INITIALIZED = False
logger = logging.getLogger("dysh")


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
    },
    "loggers": {
        "dysh": {
            "handlers": [],
            "level": "WARNING",
        },
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
