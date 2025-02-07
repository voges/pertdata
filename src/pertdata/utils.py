"""Utilities."""

import importlib.resources as pkg_resources
import json
import os
from typing import Dict

import requests
from tqdm import tqdm


def print_message(message: str, silent: bool) -> None:
    """Print a message if silent is False."""
    if not silent:
        print(message)


def download_file(url: str, file_path: str, silent: bool = False) -> None:
    """Download a file with a progress bar.

    The progress bar will display the size in binary units (e.g., KiB for kibibytes,
    MiB for mebibytes, GiB for gibibytes, etc.), which are based on powers of 1024.

    Args:
        url: The URL of the data.
        file_path: The path to save the data.
        silent: If True, do not display the progress bar and other messages.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP
            request.
        OSError: If there is an issue with writing the file.
    """
    if not os.path.exists(file_path):
        print_message(f"Downloading: {url} -> {file_path}", silent=silent)
        with requests.get(url=url, stream=True) as response:
            response.raise_for_status()
            total_size_in_bytes = int(response.headers.get("content-length", 0))
            print_message(f"Total size: {total_size_in_bytes:,} bytes", silent=silent)
            block_size = 1024
            if not silent:
                progress_bar = tqdm(
                    total=total_size_in_bytes, unit="iB", unit_scale=True
                )
            with open(file=file_path, mode="wb") as file:
                for data in response.iter_content(chunk_size=block_size):
                    if not silent:
                        progress_bar.update(len(data))
                    file.write(data)
            if not silent:
                progress_bar.close()
        print_message(f"Download completed: {file_path}", silent=silent)
    else:
        print_message(f"File already exists: {file_path}", silent=silent)


def datasets() -> Dict[str, dict]:
    """Return a dictionary of available datasets.

    The keys are the names of the datasets, and the values contain the JSON metadata.
    """
    resources_dir = pkg_resources.contents(package="pertdata.resources")
    datasets = {}
    for resource in resources_dir:
        with pkg_resources.open_text(
            package="pertdata.resources", resource=resource
        ) as json_file:
            metadata = json.load(json_file)
            name_without_extension = os.path.splitext(resource)[0]
            datasets[name_without_extension] = metadata
    return datasets
