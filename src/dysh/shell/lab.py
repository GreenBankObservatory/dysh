import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "dysh in Jupyter Lab.\n\nAll CLI arguments other than those defined below are passed through "
            "to jupyter; see $ jupyter lab --help for more details"
        )
    )

    return parser.parse_known_args()


def main():
    args, remaining_args = parse_args()

    subprocess.run(["jupyter", "lab", *remaining_args])


if __name__ == "__main__":
    main()
