"""The PertDataset class for handling perturbation datasets."""

import os
from typing import Optional

import pandas as pd
import scanpy as sc
from anndata import AnnData

import pertdata.preprocess as preprocess


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
            datasets_dir_path=dir_path,
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


def _load(dataset_name: str, dataset_variant: str, datasets_dir_path: str) -> AnnData:
    """Load perturbation dataset.

    Args:
        dataset_name: The name of the dataset.
        dataset_variant: The variant of the dataset.
        datasets_dir_path: The datasets directory path.

    Returns:
        The perturbation dataset as an AnnData object.

    Raises:
        ValueError: If the dataset name or dataset variant is unsupported.
    """
    dir_path = os.path.join(datasets_dir_path, dataset_name, dataset_variant)

    if dataset_name == "replogle_2020_exp7":
        if dataset_variant == "preprocessed":
            if not os.path.exists(path=dir_path):
                preprocess.preprocess(
                    datasets_dir_path=datasets_dir_path,
                    dataset_name=dataset_name,
                )
        else:
            raise ValueError(f"Unsupported dataset variant: {dataset_variant}")
    else:
        raise ValueError(f"Unsupported dataset name: {dataset_name}")

    # Load the dataset.
    h5ad_file_path = os.path.join(dir_path, "adata.h5ad")
    print(f"Loading: {h5ad_file_path}")
    adata = sc.read_h5ad(filename=h5ad_file_path)

    return adata
