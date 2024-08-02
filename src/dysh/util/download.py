import argparse
import sys
import uuid
from pathlib import Path

import httpx
from rich.progress import Progress


def from_url(url, path=Path(".")):
    """
    Download a file from `url` to `path`.

    Parameters
    ----------
    url : str
        The URL of the data file
    path : `pathlib.Path`
        The path to the directory to save the data. If `path/filename` already exists, the file will not be downloaded again.

    Returns
    -------
    savepath : `pathlib.Path`
        The path to the downloaded (or existing) data
    """

    try:
        # Make the HTTPX client
        client = httpx.Client(follow_redirects=True)

        with client.stream("GET", url) as resp:

            # Get the filename from the URL
            filename = Path(resp.url.path).name

            # Check if given path is directory or filename
            if path.is_dir():
                path.mkdir(parents=True, exist_ok=True)
                savepath = path / filename
            else:
                savepath = path

            # Skip downloading if file already exists
            if savepath.exists():
                print(f"{filename} already downloaded at {path}")
                return path

            else:
                # Download the file
                print(f"Attempting to download {filename}...")
                resp.raise_for_status()

                # Download to a temporary path first
                # so desired file only shows up if successful
                tmp_path = savepath.parent / f"{filename}.{uuid.uuid4()}.tmp"

                # Write chunks to file with progress bar
                with open(tmp_path, "wb") as out_file:
                    with Progress() as progress:
                        task_length = int(resp.headers.get("content-length", 0))
                        task = progress.add_task("[red]Downloading...", total=task_length)
                        for chunk in resp.iter_raw():
                            out_file.write(chunk)
                            progress.update(task, advance=len(chunk))

    # If something goes wrong, throw Exception
    except Exception as exc:
        print(exc, file=sys.stderr)
        raise

    # Rename the temp file to the desired name
    tmp_path.rename(savepath)
    print(f"Saved {filename} to {savepath}")

    return Path(savepath)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        type=Path,
        help="The path where the downloaded file will be saved. If not provided, it will be derived from the URL.",
        default=Path("."),
    )

    return parser.parse_args()


def main():
    args = parse_args()
    savepath = from_url(args.url, args.output_path)
    return savepath


if __name__ == "__main__":
    main()
