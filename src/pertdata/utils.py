"""Utility functions."""

import subprocess
from typing import Optional


def get_git_root() -> Optional[str]:
    """Return the root directory of the current Git repository.

    Returns:
        The root directory of the current Git repository, or None if the command fails.
    """
    try:
        return subprocess.check_output(
            args=["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.STDOUT,
            text=True,
        ).strip()
    except subprocess.CalledProcessError as e:
        print(f"Failed to get Git root: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return None
