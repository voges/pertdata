"""The PertDataset class for handling perturbation datasets."""

import importlib.resources as pkg_resources
import json
import os
import zipfile

import pandas as pd
import scanpy as sc
import scipy.sparse
from anndata import AnnData
from appdirs import user_cache_dir

from pertdata.utils import datasets, download_file, print_message


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
        path: The path where the dataset is stored.
        adata: The actual perturbation data.
    """

    def __init__(
        self, name: str, cache_dir_path: str = None, silent: bool = True
    ) -> "PertDataset":
        """Initialize the PertDataset object.

        Args:
            name: The name of the dataset.
            cache_dir_path: The path to the cache directory. If None, the user
                cache directory as chosen by the appdirs package is used.
            silent: If True, suppress messages.

        Returns:
            A PertDataset object.
        """
        self.name = name
        self.cache_dir_path = (
            os.path.abspath(path=cache_dir_path)
            if cache_dir_path is not None
            else user_cache_dir(appname="pertdata", appauthor=False)
        )
        self.silent = silent
        self.path = os.path.join(self.cache_dir_path, name)
        self.adata = self._load()

    def __str__(self) -> str:
        """Return a string representation of the PertDataset object."""
        return (
            f"PertDataset object\n"
            f"    name: {self.name}\n"
            f"    cache_dir_path: {self.cache_dir_path}\n"
            f"    path: {self.path}\n"
            f"    adata: AnnData object with n_obs ✕ n_vars "
            f"= {self.adata.shape[0]} ✕ {self.adata.shape[1]}"
        )

    def _load(self) -> AnnData:
        """Load perturbation dataset.

        Returns:
            The perturbation dataset as an AnnData object.
        """
        # Check if the dataset is supported.
        available_datasets = datasets()
        if f"{self.name}" not in available_datasets.keys():
            print_message("Available datasets:", silent=self.silent)
            for key in available_datasets.keys():
                print_message(f"  {key}", silent=self.silent)
            raise ValueError(f"Unsupported dataset: {self.name}")

        # Load the dataset.
        with pkg_resources.open_text(
            package="pertdata.resources", resource=f"{self.name}.json"
        ) as json_file:
            # Get the metadata.
            metadata = json.load(json_file)
            name = metadata.get("name")
            repository = metadata.get("repository")
            url = metadata.get("url")

            # Set H5AD file path.
            if repository == "GEARS":
                h5ad_file_path = os.path.join(self.path, name, "perturb_processed.h5ad")
            elif repository == "SENA" or "scPerturb":
                h5ad_file_path = os.path.join(self.path, "adata.h5ad")
            else:
                raise ValueError(f"Unsupported repository: {repository}")

            # Download the dataset if it is not already cached.
            if not os.path.exists(path=self.path):
                if repository == "GEARS":
                    os.makedirs(name=self.path, exist_ok=True)
                    zip_file_path = os.path.join(self.path, "data.zip")
                    download_file(url=url, file_path=zip_file_path, silent=self.silent)
                    with zipfile.ZipFile(file=zip_file_path, mode="r") as zip:
                        zip.extractall(path=self.path)
                elif repository == "SENA" or "scPerturb":
                    os.makedirs(name=self.path, exist_ok=True)
                    download_file(url=url, file_path=h5ad_file_path, silent=self.silent)
                else:
                    raise ValueError(f"Unsupported repository: {repository}")
            else:
                print_message(
                    f"Dataset already cached: {self.path}", silent=self.silent
                )

            # Load the dataset.
            print_message(f"Loading: {h5ad_file_path}", silent=self.silent)
            adata = sc.read_h5ad(filename=h5ad_file_path)

            return adata

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
            print_message(
                f"Exporting all {n_obs} samples to: {file_path}", silent=self.silent
            )
        elif n_samples > n_obs:
            raise ValueError(f"n_samples exceeds available samples. Max is {n_obs}.")
        else:
            print_message(
                f"Exporting the first {n_samples}/{n_obs} samples to: {file_path}",
                silent=self.silent,
            )

        # Extract cell identifiers and gene identifiers.
        cell_ids = self.adata.obs_names[:n_samples].tolist()
        gene_ids = self.adata.var_names.tolist()

        # Get the first n_samples from the expression matrix.
        expression_matrix = self.adata.X[:n_samples, :]
        if scipy.sparse.issparse(expression_matrix):
            expression_matrix = expression_matrix.todense()

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
