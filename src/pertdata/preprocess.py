"""Preprocess perturbation data."""

import argparse
import os

import scanpy as sc

import pertdata.replogle_2020 as replogle_2020


def _preprocess(datasets_dir_path: str, dataset_name: str) -> None:
    """Preprocess a dataset.

    Args:
        datasets_dir_path: The path to the datasets directory.
        dataset_name: The name of the dataset.
    """
    # Create the datasets directory. It's okay if it already exists.
    os.makedirs(name=datasets_dir_path, exist_ok=True)

    # Create the dataset directory. Raise an exception if it already exists, because
    # we don't want to overwrite an existing dataset.
    dataset_dir_path = os.path.join(datasets_dir_path, dataset_name)
    if os.path.exists(path=dataset_dir_path):
        raise FileExistsError(f"Dataset directory already exists: {dataset_dir_path}")
    os.makedirs(name=dataset_dir_path, exist_ok=True)

    # Create the "raw" directory.
    raw_dir_path = os.path.join(dataset_dir_path, "raw")
    os.makedirs(name=raw_dir_path, exist_ok=True)

    # Create the "preprocessed" directory.
    preprocessed_dir_path = os.path.join(dataset_dir_path, "preprocessed")
    os.makedirs(name=preprocessed_dir_path, exist_ok=True)

    # Download the raw data.
    if dataset_name == "adamson":
        # adamson.download_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "dixit":
        # dixit.download_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "norman":
        # norman.download_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_k562_essential":
        # replogle_k562_essential.download_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_rpe1_essential":
        # replogle_rpe1_essential.download_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_2020":
        replogle_2020.download_raw_data(dir_path=raw_dir_path)
    else:
        raise ValueError(f"Unsupported dataset: {dataset_name}")

    # Load the raw data into an AnnData object.
    print(f"Loading raw data from: {raw_dir_path}")
    if dataset_name == "adamson":
        # adata = adamson.load_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "dixit":
        # adata = dixit.load_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "norman":
        # adata = norman.load_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_k562_essential":
        # adata = replogle_k562_essential.load_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_rpe1_essential":
        # adata = replogle_rpe1_essential.load_raw_data(dir_path=raw_dir_path)
        pass
    elif dataset_name == "replogle_2020":
        adata = replogle_2020.load_raw_data(dir_path=raw_dir_path)
    else:
        raise ValueError(f"Unsupported dataset: {dataset_name}")

    # Normalize and filter the loaded data.
    sc.pp.normalize_total(adata=adata, inplace=True)
    sc.pp.log1p(adata=adata, inplace=True)
    sc.pp.highly_variable_genes(
        adata=adata, n_top_genes=5000, subset=True, inplace=True
    )

    # Save the preprocessed data to an H5AD file.
    h5ad_file_path = os.path.join(preprocessed_dir_path, "adata.h5ad")
    print(f"Saving the preprocessed data to: {h5ad_file_path}")
    adata.write(filename=h5ad_file_path)


def _print_banner(message: str) -> None:
    """Print a banner with the given message.

    Args:
        message: The message to print in the banner.
    """
    banner_width = 80
    print("\n" + "=" * banner_width)
    print(f"{message:^{banner_width}}")
    print("=" * banner_width + "\n")


def _parse_command_line_arguments() -> argparse.Namespace:
    """Parse the command line arguments.

    Returns:
        The parsed command line arguments.
    """
    dataset_names = [
        "adamson",
        "dixit",
        "norman",
        "replogle_k562_essential",
        "replogle_rpe1_essential",
        "replogle_2020",
    ]

    parser = argparse.ArgumentParser(description="Preprocess perturbation data.")

    parser.add_argument(
        "-d",
        "--datasets_dir_path",
        type=str,
        required=True,
        help="The path to the datasets directory.",
    )
    parser.add_argument(
        "-n",
        "--dataset_name",
        type=str,
        choices=dataset_names,
        required=True,
        help="The name of the dataset.",
    )

    return parser.parse_args()


def main() -> None:
    """Preprocess."""
    args = _parse_command_line_arguments()
    _print_banner(f"Preprocessing dataset: {args.dataset_name}")
    _preprocess(
        datasets_dir_path=args.datasets_dir_path,
        dataset_name=args.dataset_name,
    )


if __name__ == "__main__":
    main()
