import os
import sys
import uuid

import httpx
from rich.progress import Progress


def from_url(url, path="."):
    """
    Download a file from `url` to `path`.

    Parameters
    ----------
    url : str
        The URL of the data file
    path : str
        The path to the directory to save the data. If `path/filename` already exists, the file will not be downloaded again.

    Returns
    -------
        savepath : str
            The path to the downloaded (or existing) data
    """

    # If the URL has already been downloaded, we can skip downloading it again.
    filename = os.path.basename(url)
    savepath = f"{path}/{filename}"

    if os.path.exists(savepath):
        print(f"{filename} already downloaded at {path}")
        return path

    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    client = httpx.Client(follow_redirects=True)

    try:
        with client.stream("GET", url) as resp:
            print(f"Downloading {filename}...")
            resp.raise_for_status()

            # Download to a temporary path first
            # so desired file only shows up if successful
            tmp_path = f"{path}.{uuid.uuid4()}.tmp"

            # Write chunks to file with progress bar
            with open(tmp_path, "wb") as out_file:
                with Progress() as progress:
                    task_length = int(resp.headers.get("content-length", 0))
                    task = progress.add_task("[red]Downloading...", total=task_length)
                    for chunk in resp.iter_raw():
                        out_file.write(chunk)
                        progress.update(task, advance=len(chunk))

    # If something goes wrong, throw Exception and delete .tmp file
    except Exception as exc:
        os.remove(tmp_path)
        print(exc, file=sys.stderr)
        raise

    # Rename the temp file to the desired name
    os.rename(tmp_path, savepath)
    print(f"Saved {filename} to {savepath}")

    return savepath


if __name__ == "__main__":
    url = "http://www.gb.nrao.edu/dysh/example_data/fs-L/data/AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits"
    savepath = "/home/sandboxes/vcatlett/repos/github/GBO/dysh/"
    savepath = from_url(url, savepath)
    os.remove(savepath)
