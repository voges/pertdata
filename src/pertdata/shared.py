"""Preprocessing functionality that is shared between different datasets."""

import gzip
import os
import shutil
import tarfile
import tempfile

import anndata as ad
import pandas as pd
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


def modify_features_file(file_path: str) -> None:
    """Mofile the "<prefix>_features.tsv.gz" file.

    Modify the "<prefix>_features.tsv.gz" file to have a third column with the value
    "Gene Expression".

    Args:
        file_path: The path to the "<prefix>_features.tsv.gz" file.
    """
    with gzip.open(filename=file_path, mode="rt") as input_file:
        lines = input_file.readlines()
        if not lines[0].strip().endswith("Gene Expression"):
            with tempfile.NamedTemporaryFile() as temp_file:
                temp_file_path = temp_file.name
                with gzip.open(filename=temp_file_path, mode="wt") as output_file:
                    for line in lines:
                        line = line.strip() + "\tGene Expression\n"
                        output_file.write(line)
                shutil.copy2(src=temp_file_path, dst=file_path)


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


def filter_barcodes_and_add_condition(
    adata: ad.AnnData, barcodes_file_path: str
) -> ad.AnnData:
    """Filter the AnnData object to keep only those cells in the barcodes file.

    Filter the AnnData object to keep only cells present in the barcodes file. Also
    add the "condition" info for every cell from the barcodes file to the AnnData
    object.

    The barcodes file should have the following format:
    ```
    cell_id,condition
    AAACATACACCGAT,CREB1
    AAACATACAGAGAT,ctrl
    ...
    ```

    Args:
        adata: The AnnData object to filter.
        barcodes_file_path: The path to the barcodes file.

    Returns:
        The filtered AnnData object.
    """
    # Load the barcodes file.
    barcodes_df = pd.read_csv(filepath_or_buffer=barcodes_file_path, sep=",")

    # Get the barcodes to keep.
    barcodes_to_keep = barcodes_df["cell_id"].values

    # Filter the AnnData object.
    adata_filtered = adata[adata.obs_names.isin(values=barcodes_to_keep)].copy()

    # Add the "condition" info.
    barcode_dict = dict(zip(barcodes_df["cell_id"], barcodes_df["condition"]))
    adata_filtered.obs["condition"] = adata_filtered.obs_names.map(barcode_dict)

    return adata_filtered
