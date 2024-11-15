import argparse
import os

import pertdata.adamson as adamson
import pertdata.norman as norman
from pertdata.extract_gears_obs import extract_gears_obs
from pertdata.shared import filter_barcodes_and_add_condition


def _preprocess(datasets_dir_path: str, dataset_name: str) -> None:
    """Preprocess a dataset.

    Args:
        datasets_dir_path: The path to the datasets directory.
        dataset_name: The name of the dataset.
    """
    # Abort if the directory for the current dataset already exists. Otherwise, create
    # it.
    dataset_dir_path = os.path.join(datasets_dir_path, dataset_name)
    if os.path.exists(dataset_dir_path):
        raise FileExistsError(
            f"The dataset directory already exists: {dataset_dir_path}"
        )
    os.makedirs(name=dataset_dir_path, exist_ok=False)

    # Create the "raw" directory.
    raw_dir_path = os.path.join(datasets_dir_path, dataset_name, "raw")
    os.makedirs(name=raw_dir_path, exist_ok=True)

    # Download the raw data.
    if dataset_name == "adamson":
        adamson.download_raw_data(dir_path=raw_dir_path)
    elif dataset_name == "norman":
        norman.download_raw_data(dir_path=raw_dir_path)
    else:
        raise ValueError(f"Unsupported dataset: {dataset_name}")

    # Load the data into an AnnData object.
    print(f"Loading raw data from: {raw_dir_path}")
    if dataset_name == "adamson":
        adata = adamson.load_raw_data(dir_path=raw_dir_path)
    elif dataset_name == "norman":
        adata = norman.load_raw_data(dir_path=raw_dir_path)
    else:
        raise ValueError(f"Unsupported dataset: {dataset_name}")
    print(adata)

    apply_gears_filter = True
    if apply_gears_filter:
        # Extract the GEARS barcodes.
        gears_barcodes_file_path = extract_gears_obs(
            dataset_name=dataset_name, datasets_dir_path=datasets_dir_path
        )

        # Filter the data to keep only those cells as used by GEARS.
        print(f"Filtering the raw data based on: {gears_barcodes_file_path}")
        adata = filter_barcodes_and_add_condition(
            adata=adata, barcodes_file_path=gears_barcodes_file_path
        )
        print(adata)

    # Create the "preprocessed" directory.
    preprocessed_dir_path = os.path.join(
        datasets_dir_path, dataset_name, "preprocessed"
    )
    os.makedirs(name=preprocessed_dir_path, exist_ok=True)

    # Save the data to an H5AD file.
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
    ]

    parser = argparse.ArgumentParser(description="Preprocess the raw data.")

    parser.add_argument(
        "--datasets_dir_path",
        type=str,
        required=True,
        help="The path to the datasets directory.",
    )
    parser.add_argument(
        "--dataset_name",
        type=str,
        choices=dataset_names,
        required=True,
        help="The name of the dataset.",
    )

    return parser.parse_args()


def main() -> None:
    """Preprocess the raw data."""
    args = _parse_command_line_arguments()
    _print_banner(f"Preprocessing dataset: {args.dataset_name}")
    _preprocess(
        datasets_dir_path=args.datasets_dir_path,
        dataset_name=args.dataset_name,
    )


if __name__ == "__main__":
    main()
