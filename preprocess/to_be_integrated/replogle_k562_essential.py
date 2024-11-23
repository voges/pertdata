"""Preprocessing for the Replogle K562 Essential dataset."""

import os
import sys

import scanpy as sc

# Add the root of the project to sys.path.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from preprocess.shared import download_file


def download_raw_data(dir_path: str) -> None:
    """Download the raw data.

    Args:
        dir_path: The directory path where the raw data will be stored.
    """
    download_file(
        url="https://plus.figshare.com/ndownloader/files/35773219",
        path=os.path.join(dir_path, "K562_essential_raw_singlecell_01.h5ad"),
    )


def load_raw_data(dir_path: str) -> sc.AnnData:
    """Load the raw data.

    Args:
        dir_path: The directory path where the raw data is stored.

    Returns:
        The AnnData object.
    """
    adata = sc.read_h5ad(
        filename=os.path.join(dir_path, "K562_essential_raw_singlecell_01.h5ad")
    )
    return adata
