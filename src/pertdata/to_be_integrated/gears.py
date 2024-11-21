"""Utility to extract ['cell_id', 'condition'] from the GEARS datasets."""

import os
import sys

import pandas as pd

# Add the root of the project to sys.path.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import transmet.pert_dataset as pt


def _fix_adamson(csv_file_path: str) -> None:
    """Fix the CSV file for the Adamson dataset.

    Args:
        csv_file_path: The path to the CSV file.
    """
    # Read the CSV file into a DataFrame.
    df = pd.read_csv(filepath_or_buffer=csv_file_path)

    # Remove the "(?)" suffix from all rows.
    df = df.replace(to_replace=r"\(\?\)", value="", regex=True)

    # Save the modified DataFrame back to the same CSV file.
    df.to_csv(path_or_buf=csv_file_path, index=False)


def extract_gears_obs(dataset_name: str, datasets_dir_path: str) -> str:
    """Extract ['cell_id', 'condition'] from the GEARS datasets.

    Args:
        dataset_name: The name of the dataset.
        datasets_dir_path: The path to the datasets directory.

    Returns:
        The path to the CSV file containing the observations.
    """
    # Make the output file path.
    output_file_path = os.path.join(datasets_dir_path, dataset_name, "gears", "obs.csv")

    if not os.path.exists(output_file_path):
        # Load the dataset.
        pert_dataset = pt.PertDataset(
            name=dataset_name, variant="gears", dir_path=datasets_dir_path
        )

        # Export the observations.
        print(f"Exporting observations to: {output_file_path}")
        pert_dataset.export_obs_to_csv(
            file_path=output_file_path, obs=["cell_id", "condition"]
        )

        # Apply fixes.
        if dataset_name == "adamson":
            _fix_adamson(csv_file_path=output_file_path)

    return output_file_path


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


def apply_gears_filter():
    # Extract the GEARS barcodes.
    gears_barcodes_file_path = extract_gears_obs(
        dataset_name=dataset_name, datasets_dir_path=datasets_dir_path
    )

    # Filter the data to keep only those cells as used by GEARS.
    print(f"Filtering the raw data based on: {gears_barcodes_file_path}")
    adata = filter_barcodes_and_add_condition(
        adata=adata, barcodes_file_path=gears_barcodes_file_path
    )


"""Utility to extract ['cell_id', 'condition'] from the GEARS datasets."""

import gzip
import os
import shutil
from typing import List

import gdown
import scanpy as sc
from anndata import AnnData


def _fix_adamson(csv_file_path: str) -> None:
    """Fix the CSV file for the Adamson dataset.

    Args:
        csv_file_path: The path to the CSV file.
    """
    # Read the CSV file into a DataFrame.
    df = pd.read_csv(filepath_or_buffer=csv_file_path)

    # Remove the "(?)" suffix from all rows.
    df = df.replace(to_replace=r"\(\?\)", value="", regex=True)

    # Save the modified DataFrame back to the same CSV file.
    df.to_csv(path_or_buf=csv_file_path, index=False)


def _export_obs_to_csv(adata: AnnData, file_path: str, obs: List[str]) -> None:
    """Export the observation data to a CSV file.

    Args:
        adata: The AnnData object.
        file_path: The path to save the CSV file.
        obs: The list of observations to export.
    """
    # Get a copy of the observation data.
    all_obs = adata.obs.copy()

    # Handle special case: "cell_id" is the obs index. Include it as the first
    # column of our temporary DataFrame to be able to export it.
    if "cell_id" in obs:
        all_obs["cell_id"] = adata.obs.index
        all_obs = all_obs[
            ["cell_id"] + [item for item in obs if item != "cell_id"]
        ]  # Reorder columns.

    # Check if the requested observations are present in the dataset.
    for o in obs:
        if o not in all_obs.columns:
            raise ValueError(f"Observation '{o}' not found in the dataset.")

    # Make a DataFrame with only the requested observations.
    requested_obs = all_obs[obs]

    # Export the observation data to a CSV file.
    requested_obs[obs].to_csv(path_or_buf=file_path, index=False)


def extract_gears_obs(dataset_name: str, datasets_dir_path: str) -> str:
    """Extract ['cell_id', 'condition'] from a GEARS dataset.

    Args:
        dataset_name: The name of the dataset.
        datasets_dir_path: The path to the datasets directory.

    Returns:
        The path to the CSV file containing the observations.
    """
    # Make the output file path.
    output_file_path = os.path.join(datasets_dir_path, dataset_name, "gears", "obs.csv")

    if not os.path.exists(output_file_path):
        if dataset_name == "adamson":
            url = "https://drive.google.com/uc?id=1W1phErDoQ9U5iJZSEuyEZM4dF8U8ZWqf"
        elif dataset_name == "dixit":
            url = "https://drive.google.com/uc?id=1BN6gwKFgJIpR9fXfdQ9QeHm8mAzvmhKQ"
        elif dataset_name == "norman":
            url = "https://drive.google.com/uc?id=1T5_varFOGWUtSig4RQRCSsfPxivUwd9j"
        elif dataset_name == "replogle_k562_essential":
            url = "https://drive.google.com/uc?id=12flxmpj-XnJ8BZKtf-sgBhdUN2X4v7CD"
        elif dataset_name == "replogle_rpe1_essential":
            url = "https://drive.google.com/uc?id=1b-ZwE_Y6dNKqb4KQgUgFKfl6OGC8lmYE"
        else:
            raise ValueError(f"Unsupported dataset name: {dataset_name}")

        # Make the directory for the selected dataset.
        dataset_dir_path = os.path.join(datasets_dir_path, dataset_name, "gears")
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

        # Export the observations.
        print(f"Exporting observations to: {output_file_path}")
        _export_obs_to_csv(
            adata=adata, file_path=output_file_path, obs=["cell_id", "condition"]
        )

        # Apply fixes.
        if dataset_name == "adamson":
            _fix_adamson(csv_file_path=output_file_path)

    return output_file_path
