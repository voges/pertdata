"""Functionality for handling perturbation datasets."""

import gzip
import os
import shutil
from typing import Optional

import gdown
import pandas as pd
import scanpy as sc
from anndata import AnnData


class PertDataset:
    """Class for handling a perturbation dataset.

    The actual perturbation dataset is stored in an
    [AnnData](https://anndata.readthedocs.io/en/latest/) object.

    AnnData is specifically designed for matrix-like data. By this we mean that we have
    n observations, each of which can be represented as d-dimensional vectors, where
    each dimension corresponds to a variable or feature. Both the rows and columns of
    this matrix are special in the sense that they are indexed.

    For instance, in scRNA-seq data:
    - Each row corresponds to a cell with a cell identifier (i.e., barcode).
    - Each column corresponds to a gene with a gene identifier.

    Attributes:
        name: The name of the dataset.
        variant: The variant of the dataset.
        path: The path where the dataset is stored.
        adata: The actual perturbation data.
    """

    def __init__(self, name: str, variant: str, dir_path: str) -> "PertDataset":
        """Initialize the PertDataset object.

        Args:
            name: The name of the dataset.
            variant: The variant of the dataset.
            dir_path: The path to the datasets directory.

        Returns:
            A PertDataset object.
        """
        # Initialize the attributes.
        self.name: Optional[str] = None
        self.variant: Optional[str] = None
        self.path: Optional[str] = None
        self.adata: Optional[AnnData] = None

        # Set the attributes.
        self.name = name
        self.variant = variant
        self.path = os.path.join(dir_path, name, variant)
        self.adata = _load(
            dataset_name=name,
            dataset_variant=variant,
            dataset_dir_path=self.path,
        )

    def __str__(self) -> str:
        """Return a string representation of the PertDataset object."""
        return (
            f"PertDataset object\n"
            f"    name: {self.name}\n"
            f"    variant: {self.variant}\n"
            f"    path: {self.path}\n"
            f"    adata: AnnData object with n_obs ✕ n_vars "
            f"= {self.adata.shape[0]} ✕ {self.adata.shape[1]}"
        )

    def normalize_(self, type: str = "CPM") -> None:
        """Normalize the gene expression matrix.

        Args:
            type: The type of normalization to apply (supported: "CPM").
        """
        if type == "CPM":
            _normalize_cpm_(adata=self.adata)
        else:
            raise ValueError(f"Unsupported normalization type: {type}")

    def export_tsv(self, file_path: str, n_samples: int = None) -> None:
        """Save the perturbation data to a TSV file.

        If n_samples is provided, only the first n_samples samples are exported.

        The TSV file has the following format:
        - The first row contains the cell identifiers.
        - The first column contains the gene identifiers.
        - The remaining entries are the gene expression values.

        Args:
            file_path: The path to save the TSV file.
            n_samples: The number of samples to export.
        """
        # Export all samples if n_samples is not provided.
        n_obs = self.adata.shape[0]
        if n_samples is None:
            n_samples = n_obs
            print(f"Exporting all {n_obs} samples to: {file_path}")
        elif n_samples > n_obs:
            raise ValueError(f"n_samples exceeds available samples. Max is {n_obs}.")
        else:
            print(f"Exporting the first {n_samples}/{n_obs} samples to: {file_path}")

        # Extract cell identifiers and gene identifiers.
        cell_ids = self.adata.obs_names[:n_samples].tolist()
        gene_ids = self.adata.var_names.tolist()

        # Get the first n_samples from the expression matrix.
        expression_matrix = self.adata.X[:n_samples, :].todense()

        # Transpose expression matrix to match the desired output (genes as rows,
        # cells as columns).
        expression_matrix = expression_matrix.T

        # Create a DataFrame for export.
        expression_df = pd.DataFrame(
            data=expression_matrix, index=gene_ids, columns=cell_ids
        )

        # Reset index to move the gene identifiers (row index) to a column.
        expression_df.reset_index(inplace=True)

        # Rename the index column to "Gene" for clarity.
        expression_df.rename(columns={"index": "Gene"}, inplace=True)

        # Save the DataFrame to a TSV file.
        expression_df.to_csv(path_or_buf=file_path, sep="\t", index=False)


def _load(dataset_name: str, dataset_variant: str, dataset_dir_path: str) -> AnnData:
    """Load perturbation dataset.

    Args:
        dataset_name: The name of the dataset.
        dataset_variant: The variant of the dataset.
        dataset_dir_path: The directory path where the dataset is stored.

    Returns:
        The perturbation dataset as an AnnData object.

    Raises:
        ValueError: If the dataset name or dataset variant is unsupported.
    """
    if dataset_name == "adamson":
        if dataset_variant == "gears":
            url = "https://drive.google.com/uc?id=1W1phErDoQ9U5iJZSEuyEZM4dF8U8ZWqf"
        elif dataset_variant == "preprocessed":
            url = "https://drive.google.com/uc?id=150AH0uEJYDERCHnnM2P9knZpmY4tOGNZ"
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    elif dataset_name == "dixit":
        if dataset_variant == "gears":
            url = "https://drive.google.com/uc?id=1BN6gwKFgJIpR9fXfdQ9QeHm8mAzvmhKQ"
        elif dataset_variant == "preprocessed":
            url = "https://drive.google.com/uc?id=1qVQD9zU_Dpj5AUJaD3tM8BpfTkKqgvog"
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    elif dataset_name == "norman":
        if dataset_variant == "gears":
            url = "https://drive.google.com/uc?id=1T5_varFOGWUtSig4RQRCSsfPxivUwd9j"
        elif dataset_variant == "preprocessed":
            url = "https://drive.google.com/uc?id=1KgatnkMHy-3uh5tNkAFrsmGrivDOUjpx"
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    elif dataset_name == "replogle_k562_essential":
        if dataset_variant == "gears":
            url = "https://drive.google.com/uc?id=12flxmpj-XnJ8BZKtf-sgBhdUN2X4v7CD"
        elif dataset_variant == "preprocessed":
            url = "https://drive.google.com/uc?id=1zdznT5x92pd8-o9HOtudKkwvj4CgqnNO"
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    elif dataset_name == "replogle_rpe1_essential":
        if dataset_variant == "gears":
            url = "https://drive.google.com/uc?id=1b-ZwE_Y6dNKqb4KQgUgFKfl6OGC8lmYE"
        elif dataset_variant == "preprocessed":
            url = "https://drive.google.com/uc?id=1KhmS61iYwE0ineaiQyqX3chjJmCfmtFE"
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    else:
        raise ValueError(f"Unsupported dataset name: {dataset_name}")

    # Make the directory for the selected dataset.
    os.makedirs(name=dataset_dir_path, exist_ok=True)

    # Make the file path for the selected dataset.
    h5ad_file_path = os.path.join(dataset_dir_path, "adata.h5ad")

    # Download and extract the dataset if it does not exist.
    if not os.path.exists(path=h5ad_file_path):
        # Make the file path for the compressed dataset.
        zip_file_path = os.path.join(dataset_dir_path, "adata.h5ad.gz")

        # Download the dataset.
        print(f"Downloading: {url} -> {zip_file_path}")
        gdown.download(url=url, output=zip_file_path)

        # Extract the dataset.
        print(f"Extracting: {zip_file_path} -> {h5ad_file_path}")
        with gzip.open(filename=zip_file_path, mode="rb") as f_in:
            with open(file=h5ad_file_path, mode="wb") as f_out:
                shutil.copyfileobj(fsrc=f_in, fdst=f_out)

        # Remove the compressed file.
        print(f"Removing: {zip_file_path}")
        os.remove(path=zip_file_path)

    # Load the dataset.
    print(f"Loading: {h5ad_file_path}")
    adata = sc.read_h5ad(filename=h5ad_file_path)

    return adata


def _normalize_cpm_(adata: AnnData) -> None:
    """Normalize the data (target_sum=1e6 is CPM normalization).

    Args:
        adata: The AnnData object to normalize.
    """
    sc.pp.normalize_total(adata=adata, target_sum=1e6, inplace=True)
