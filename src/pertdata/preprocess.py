"""Preprocess perturbation data."""

import argparse
import os

import pertdata.replogle_2020_exp7 as replogle_2020_exp7


def preprocess(datasets_dir_path: str, dataset_name: str) -> None:
    """Preprocess a dataset.

    Args:
        datasets_dir_path: The path to the datasets directory.
        dataset_name: The name of the dataset.
    """
    # Create the datasets directory. It's okay if it already exists.
    os.makedirs(name=datasets_dir_path, exist_ok=True)

    # Create the dataset directory. It's okay if it already exists.
    dataset_dir_path = os.path.join(datasets_dir_path, dataset_name)
    os.makedirs(name=dataset_dir_path, exist_ok=True)

    # Create the "raw" directory and download the raw data. If the "raw" directory
    # already exists, skip the download.
    raw_dir_path = os.path.join(dataset_dir_path, "raw")
    if not os.path.exists(path=raw_dir_path):
        os.makedirs(name=raw_dir_path, exist_ok=True)
        if dataset_name == "norman":
            # norman.download_raw_data(dir_path=raw_dir_path)
            pass
        elif dataset_name == "replogle_2020_exp7":
            replogle_2020_exp7.download_raw_data(dir_path=raw_dir_path)
        else:
            raise ValueError(f"Unsupported dataset: {dataset_name}")
    else:
        print(f"Skipping downloading the raw data to: {raw_dir_path}")

    # If the "preprocessed" directory does not already exists:
    # - Load the raw data into an AnnData object.
    # - Save the preprocessed data.
    preprocessed_dir_path = os.path.join(dataset_dir_path, "preprocessed")
    if not os.path.exists(path=preprocessed_dir_path):
        print(f"Loading raw data from: {raw_dir_path}")
        if dataset_name == "norman":
            # adata = norman.load_raw_data(dir_path=raw_dir_path)
            pass
        elif dataset_name == "replogle_2020_exp7":
            adata = replogle_2020_exp7.load_raw_data(dir_path=raw_dir_path)
        else:
            raise ValueError(f"Unsupported dataset: {dataset_name}")

        os.makedirs(name=preprocessed_dir_path)
        h5ad_file_path = os.path.join(preprocessed_dir_path, "adata.h5ad")
        print(f"Saving the preprocessed data to: {h5ad_file_path}")
        adata.write(filename=h5ad_file_path)
    else:
        print(f"Skipping saving the preprocessed data to: {preprocessed_dir_path}")


def _parse_command_line_arguments() -> argparse.Namespace:
    """Parse the command line arguments.

    Returns:
        The parsed command line arguments.
    """
    dataset_names = [
        # "adamson",
        # "dixit",
        "norman",
        # "replogle_k562_essential",
        # "replogle_rpe1_essential",
        "replogle_2020_exp7",
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
    preprocess(
        datasets_dir_path=args.datasets_dir_path,
        dataset_name=args.dataset_name,
    )


if __name__ == "__main__":
    main()
