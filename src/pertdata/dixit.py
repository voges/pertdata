"""Preprocessing for the Dixit dataset."""

import os
import sys

import anndata as ad
import pandas as pd
import scanpy as sc

# Add the root of the project to sys.path.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from preprocess.shared import download_and_extract_tar_file


def download_raw_data(dir_path: str) -> None:
    """Download the raw data.

    Args:
        dir_path: The directory path where the raw data will be stored.
    """
    download_and_extract_tar_file(
        url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90063&format=file",
        dir_path=dir_path,
    )


def load_raw_data(dir_path: str) -> sc.AnnData:
    """Load the raw data.

    There are several experiments in the Dixit dataset. GEARS only uses two of the
    experiments, whose files have the following prefixes:
        1. GSM2396858_k562_tfs_7
        2. GSM2396861_k562_ccycle

    Args:
        dir_path: The directory path where the raw data is stored.

    Returns:
        The AnnData object.
    """
    adata_tfs_7 = _load_raw_data_experiment(
        prefix="GSM2396858_k562_tfs_7", dir_path=dir_path
    )
    adata_ccycle = _load_raw_data_experiment(
        prefix="GSM2396861_k562_ccycle", dir_path=dir_path
    )
    adata_combined = ad.concat([adata_tfs_7, adata_ccycle], join="outer")
    return adata_combined


def _load_raw_data_experiment(prefix: str, dir_path: str) -> sc.AnnData:
    """Load the raw data for a specific experiment.

    Args:
        prefix: The prefix of the files for the experiment.
        dir_path: The directory path where the raw data is stored.

    Returns:
        The AnnData object.
    """
    # Make the file paths.
    mtx_file_path = os.path.join(dir_path, f"{prefix}.mtx.txt.gz")
    genenames_file_path = os.path.join(dir_path, f"{prefix}_genenames.csv.gz")
    cellnames_file_path = os.path.join(dir_path, f"{prefix}_cellnames.csv.gz")

    # Load the data.
    adata = sc.read_mtx(filename=mtx_file_path).transpose()

    # Set the gene identifiers.
    gene_names = pd.read_csv(filepath_or_buffer=genenames_file_path, index_col=0)
    gene_ids = gene_names["0"].str.split("_").str[0]  # Extract only the gene IDs.
    adata.var_names = gene_ids

    # Set the cell identifiers.
    cell_names = pd.read_csv(filepath_or_buffer=cellnames_file_path, index_col=0)
    cell_ids = cell_names["0"]
    adata.obs_names = cell_ids

    return adata
