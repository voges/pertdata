"""Utilities."""

import importlib.resources as pkg_resources
import json
import os
from typing import Dict

import requests
import toml
from appdirs import user_cache_dir
from tqdm import tqdm


def download_file(url: str, path: str, skip_if_exists: bool = True) -> None:
    """Download a file with a progress bar.

    The progress bar will display the size in binary units (e.g., KiB for kibibytes,
    MiB for mebibytes, GiB for gibibytes, etc.), which are based on powers of 1024.

    Args:
        url: The URL of the file.
        path: The path where the file will be saved.
        skip_if_exists: If True, skip downloading the file if it already exists.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP
            request.
        OSError: If there is an issue with writing the file.
    """
    if skip_if_exists and os.path.exists(path=path):
        print(f"Skipping download because file already exists: {path}")
        return

    print(f"Downloading: {url} -> {path}")
    progress_bar = None
    try:
        with requests.get(url=url, stream=True) as response:
            response.raise_for_status()
            total_size_in_bytes = int(
                response.headers.get(key="content-length", default=0)
            )
            print(f"Total size: {total_size_in_bytes:,} bytes")
            block_size = 1024
            with tqdm(
                total=total_size_in_bytes, unit="iB", unit_scale=True
            ) as progress_bar:
                with open(file=path, mode="wb") as file:
                    for data in response.iter_content(chunk_size=block_size):
                        progress_bar.update(n=len(data))
                        file.write(data)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        raise
    except OSError as e:
        print(f"Error writing file: {e}")
        raise


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


def cache_dir_path() -> str:
    """Return the path to the cache directory."""
    return user_cache_dir(appname="pertdata", appauthor=False)


def get_version() -> str:
    """Get the version from pyproject.toml."""
    pyproject_file_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "pyproject.toml"
    )
    with open(file=pyproject_file_path, mode="r") as pyproject_file:
        pyproject_data = toml.load(pyproject_file)
        return pyproject_data.get("project", {}).get("version")
