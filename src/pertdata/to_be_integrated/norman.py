"""Preprocessing for the Norman dataset."""

import gzip
import os
import shutil
import tempfile

import scanpy as sc

from pertdata.shared import download_file


def _modify_features_file(file_path: str) -> None:
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


def download_raw_data(dir_path: str) -> None:
    """Download the raw data.

    Args:
        dir_path: The directory path where the raw data will be stored.
    """
    urls = [
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133344&format=file&file=GSE133344%5Fraw%5Fbarcodes%2Etsv%2Egz",
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133344&format=file&file=GSE133344%5Fraw%5Fcell%5Fidentities%2Ecsv%2Egz",
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133344&format=file&file=GSE133344%5Fraw%5Fgenes%2Etsv%2Egz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133344/suppl/GSE133344%5Fraw%5Fmatrix%2Emtx%2Egz",
    ]

    filenames = [
        "raw_barcodes.tsv.gz",
        "raw_cell_identities.csv.gz",
        "raw_genes.tsv.gz",
        "raw_matrix.mtx.gz",
    ]

    for url, filename in zip(urls, filenames):
        download_file(url=url, path=os.path.join(dir_path, filename))


def load_raw_data(dir_path: str) -> sc.AnnData:
    """Load the raw data.

    Args:
        dir_path: The directory path where the raw data is stored.

    Returns:
        The AnnData object.
    """
    # We need to rename "raw_genes.tsv.gz" to "raw_features.tsv.gz", because the
    # function sc.read_10x_mtx() expects the file to be named "raw_features.tsv.gz".
    # We make a copy and keep "raw_genes.tsv.gz" to avoid duplicate downloads.
    shutil.copy2(
        src=os.path.join(dir_path, "raw_genes.tsv.gz"),
        dst=os.path.join(dir_path, "raw_features.tsv.gz"),
    )

    # Also, the "raw_features.tsv.gz" file needs to have a third column with the value
    # "Gene Expression".
    raw_features_file_path = os.path.join(dir_path, "raw_features.tsv.gz")

    _modify_features_file(file_path=raw_features_file_path)

    adata = sc.read_10x_mtx(
        path=dir_path, var_names="gene_ids", cache=False, prefix="raw_"
    )

    return adata
