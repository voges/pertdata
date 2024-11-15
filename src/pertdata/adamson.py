"""Preprocessing for the Adamson dataset."""

import os
import shutil

import anndata as ad
import scanpy as sc

from pertdata.shared import download_and_extract_tar_file, modify_features_file


def download_raw_data(dir_path: str) -> None:
    """Download the raw data.

    Args:
        dir_path: The directory path where the raw data will be stored.
    """
    download_and_extract_tar_file(
        url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90546&format=file",
        dir_path=dir_path,
    )

    # The raw data of the Adamson dataset needs to be organized quite specifically to
    # be loaded correctly by Scanpy.
    _organize_raw_data_adamson(dir_path=dir_path)

    # The "raw_features.tsv.gz" files needs to have a third column with the value
    # "Gene Expression".
    for dir_name in os.listdir(path=dir_path):
        curr_dir_path = os.path.join(dir_path, dir_name)
        if os.path.isdir(curr_dir_path):
            raw_features_file_path = os.path.join(curr_dir_path, "raw_features.tsv.gz")
            modify_features_file(file_path=raw_features_file_path)


def load_raw_data(dir_path: str) -> sc.AnnData:
    """Load the raw data.

    Args:
        dir_path: The directory path where the raw data is stored.

    Returns:
        The AnnData object.
    """
    # Read the datasets.
    adata1 = sc.read_10x_mtx(
        path=os.path.join(dir_path, "GSM2406675_10X001"),
        var_names="gene_ids",
        cache=False,
        prefix="raw_",
    )
    adata2 = sc.read_10x_mtx(
        path=os.path.join(dir_path, "GSM2406677_10X005"),
        var_names="gene_ids",
        cache=False,
        prefix="raw_",
    )
    adata3 = sc.read_10x_mtx(
        path=os.path.join(dir_path, "GSM2406681_10X010"),
        var_names="gene_ids",
        cache=False,
        prefix="raw_",
    )

    # Concatenate the datasets with unique observations.
    adata_combined = ad.concat(
        [adata1, adata2, adata3], axis=0, join="outer", index_unique=None
    )

    return adata_combined


def _organize_raw_data_adamson(dir_path: str) -> None:
    """Organize the raw data for the Adamson dataset.

    Organize the raw data so that it can be loaded correctly by Scanpy. We need to
    create a separate directory for each experiment and move the files into the
    appropriate directories.

    Args:
        dir_path: The directory path where the raw data is stored.
    """
    # List all files in the raw data directory.
    files = os.listdir(path=dir_path)

    # Iterate through the files to identify the experiments.
    experiments = set()
    for file in files:
        if file.startswith("GSM"):
            # Extract the experiment name.
            experiment = "_".join(file.split("_")[:2])
            experiments.add(experiment)

    # For each experiment, create a directory and move the respective files into it.
    for experiment in experiments:
        # Create a new directory for the experiment.
        experiment_dir_path = os.path.join(dir_path, experiment)
        os.makedirs(name=experiment_dir_path, exist_ok=True)

        # Define the new file names.
        new_file_names = {
            "barcodes": "raw_barcodes.tsv.gz",
            "matrix": "raw_matrix.mtx.gz",
            "genes": "raw_features.tsv.gz",
            "cell_identities": "raw_cell_identities.csv.gz",
        }

        # Move (and rename) the files into the new directory.
        for file_type, new_file_name in new_file_names.items():
            for file in files:
                if file.startswith(experiment) and file_type in file:
                    old_file_path = os.path.join(dir_path, file)
                    new_file_path = os.path.join(experiment_dir_path, new_file_name)
                    shutil.move(src=old_file_path, dst=new_file_path)
