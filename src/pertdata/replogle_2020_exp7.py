"""Preprocessing for the "Replogle 2020 Exp. 7" dataset."""

# Paper: https://doi.org/10.1038/s41587-020-0470-y
# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194

import gzip
import os
import tarfile

import pandas as pd
import scanpy as sc

from pertdata.shared import download_file


def download_raw_data(dir_path: str) -> None:  # noqa: D103
    download_file(
        url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146194&format=file",
        path=os.path.join(dir_path, "GSE146194_RAW.tar"),
    )


def load_raw_data(dir_path: str) -> None:  # noqa: D103
    # # Extract the tar file.
    tar_file_path = os.path.join(dir_path, "GSE146194_RAW.tar")
    extract_dir_path = os.path.join(dir_path, "GSE146194_RAW")
    with tarfile.open(name=tar_file_path, mode="r") as tar:
        tar.extractall(path=extract_dir_path)

    # Load the data.
    adata = sc.read_10x_mtx(
        path=extract_dir_path,
        var_names="gene_ids",
        cache=False,
        prefix="GSM4367985_exp7.",
    )

    # Add info from the "*_cell_identities.csv.gz" file to the AnnData object.
    cell_identities_file_path = os.path.join(
        extract_dir_path, "GSM4367985_exp7.cell_identities.csv.gz"
    )
    barcodes_file_path = os.path.join(
        extract_dir_path, "GSM4367985_exp7.barcodes.tsv.gz"
    )
    with gzip.open(cell_identities_file_path, mode="r") as cell_identities_file:
        cell_identities_df = pd.read_csv(filepath_or_buffer=cell_identities_file)
        with gzip.open(barcodes_file_path, mode="r") as barcodes_file:
            barcodes_df = pd.read_csv(
                filepath_or_buffer=barcodes_file, header=None, names=["cell_barcode"]
            )
            merged_df = pd.merge(
                left=barcodes_df,
                right=cell_identities_df,
                on="cell_barcode",
                how="left",
            )

    # Ensure the merged_df index matches the obs_names of adata.
    merged_df.set_index("cell_barcode", inplace=True)

    # Convert all columns to strings.
    merged_df = merged_df.astype(str)

    # Add the merged_df as obs to adata.
    adata.obs = merged_df

    return adata
