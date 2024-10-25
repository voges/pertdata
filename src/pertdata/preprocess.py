import argparse
import os

from pertdata.utils import get_git_root


def _preprocess(datasets_dir_path: str, dataset_name: str) -> None:
    """Preprocess a dataset.

    Args:
        datasets_dir_path: The path to the datasets directory.
        dataset_name: The name of the dataset.
        apply_gears_filter: Whether to reduce the data to the same set of cells as used
            by GEARS.
    """
    # Create the "raw" directory.
    raw_dir_path = os.path.join(datasets_dir_path, dataset_name, "raw")
    os.makedirs(name=raw_dir_path, exist_ok=True)


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
        default=os.path.join(get_git_root(), "datasets"),
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
