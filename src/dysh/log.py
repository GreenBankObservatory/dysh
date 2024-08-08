import logging
import logging.config
import os
from datetime import datetime
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


config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "verbose": {
            "format": "%(asctime)s - %(module)s.%(funcName)s - %(levelname)s - %(message)s",
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
    },
    "loggers": {
        "dysh": {
            "handlers": ["stderr", "dysh_global_log_file"],
            "level": "WARNING",
        },
    },
}


def init_logging():
    global logger
    logging.config.dictConfig(config)
    logger = logging.getLogger("dysh")


def config_logging(verbosity: int, level: Union[int | None] = None, path: Union[Path | None] = None):
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

    logging.getLogger().setLevel(level)
    logger.setLevel(level)
    logger.debug(f"Logging has been set to verbosity {verbosity} / level {logging.getLevelName(level)}")

    if path:
        enable_instance_log(level=level, path=path)
        logger.info(f"Log file for this instance of dysh: {path}")


def enable_instance_log(level: int, path: Path):
    handler = logging.FileHandler(path)
    handler.setLevel(level)
    logger.addHandler(handler)


logger = None
