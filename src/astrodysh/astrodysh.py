"""astrodysh entry point"""

import argparse
import logging
from pathlib import Path

import pyspeckit


logger = logging.getLogger(__name__)


def init_logging(verbosity):
    if verbosity >= 2:
        level = logging.DEBUG
    elif verbosity == 1:
        level = logging.INFO
    elif verbosity == 0:
        level = logging.WARNING
    else:
        raise ValueError(f"Invalid verbosity level: {verbosity}")

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    logger.setLevel(level)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    # See: https://docs.python.org/3/library/logging.html#logrecord-attributes
    formatter = logging.Formatter("[%(asctime)s - %(levelname)s] %(message)s")
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)


def astrodysh(input_path: Path, output_path: Path):
    logger.info("hello world")
    pyspeckit.__file__
    return True


def main():
    args = parse_args()
    init_logging(args.verbosity)
    astrodysh(input_path=args.input_path, output_path=args.output_path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path, default=Path("."))
    parser.add_argument("-v", "--verbosity", type=int, choices=[0, 1, 2, 3], default=1)

    return parser.parse_args()


if __name__ == "__main__":
    main()
