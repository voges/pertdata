"""Preprocessing for the "Replogle 2020" dataset."""

# Paper: https://doi.org/10.1038/s41587-020-0470-y
# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194

import os
import tarfile

import scanpy as sc

from pertdata.shared import download_file


def download_raw_data(dir_path: str) -> None:
    """Download the raw data.

    Args:
        dir_path: The directory path where the raw data will be stored.
    """
    download_file(
        url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146194&format=file",
        path=os.path.join(dir_path, "GSE146194_RAW.tar"),
    )


def load_raw_data(dir_path: str) -> None:
    """Load the raw data into an AnnData object.

    Args:
        dir_path: The directory path where the raw data is stored.
    """
    # Extract the tar file.
    tar_file_path = os.path.join(dir_path, "GSE146194_RAW.tar")
    extract_dir_path = os.path.join(dir_path, "GSE146194_RAW")
    with tarfile.open(name=tar_file_path, mode="r") as tar:
        tar.extractall(path=extract_dir_path)

    # Fix: Rename all files with the extension ".matrix.mtx.gz" to ".mtx.gz".
    for root, _, files in os.walk(top=extract_dir_path):
        for file in files:
            if file.endswith(".matrix.mtx.gz"):
                os.rename(
                    os.path.join(root, file),
                    os.path.join(root, file.replace(".matrix.mtx.gz", ".mtx.gz")),
                )

    # Load the data.
    adata = sc.read_10x_mtx(
        path=extract_dir_path,
        var_names="gene_ids",
        cache=False,
        prefix="GSM4367989_exp11.",
    )

    return adata
