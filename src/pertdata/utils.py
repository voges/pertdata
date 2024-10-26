"""Utility functions."""

import os
import subprocess
import tarfile
from typing import Optional

import requests
from tqdm import tqdm


def download_file(url: str, path: str, skip_if_exists=True) -> None:
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
    if skip_if_exists and os.path.exists(path):
        return

    print(f"Downloading: {url} -> {path}")
    try:
        with requests.get(url=url, stream=True) as response:
            response.raise_for_status()
            total_size_in_bytes = int(
                response.headers.get(key="content-length", default=0)
            )
            print(f"Total size: {total_size_in_bytes:,} bytes")
            block_size = 1024
            progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
            with open(file=path, mode="wb") as file:
                for data in response.iter_content(chunk_size=block_size):
                    progress_bar.update(n=len(data))
                    file.write(data)
            progress_bar.close()
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
    except OSError as e:
        print(f"Error writing file: {e}")
    finally:
        if "progress_bar" in locals():
            progress_bar.close()


def download_and_extract_tar_file(url: str, dir_path: str) -> None:
    """Download and extract a TAR file.

    Args:
        url: The URL of the TAR file.
        dir_path: The directory path where the TAR file will be stored and extracted.
    """
    # Download the TAR file.
    tar_file_path = os.path.join(dir_path, "data.tar")
    download_file(url=url, path=tar_file_path)

    # Extract the TAR file.
    with tarfile.open(name=tar_file_path, mode="r") as tar_file:
        tar_file.extractall(path=dir_path, filter=None)


def get_git_root() -> Optional[str]:
    """Return the root directory of the current Git repository.

    Returns:
        The root directory of the current Git repository, or None if the command fails.
    """
    try:
        return subprocess.check_output(
            args=["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.STDOUT,
            text=True,
        ).strip()
    except subprocess.CalledProcessError as e:
        print(f"Failed to get Git root: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return None
