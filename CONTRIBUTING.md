# Contributing to `pertdata`

## Development

To install the `pertdata` package for development, follow these steps:

1. Create a [virtual environment](https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-and-using-virtual-environments).

2. Install the required dependencies:
    ```shell
    pip3 install --requirement requirements.txt
    ```

3. Install the `pertdata` package in editable mode:
    ```shell
    pip3 install --editable .
    ```

## Distribution

1. Update the version in [pyproject.toml](pyproject.toml).
    Use [Semantic Versioning](https://semver.org).
    Given a version number MAJOR.MINOR.PATCH, increment the:
    - MAJOR version when you make incompatible API changes,
    - MINOR version when you add functionality in a backward compatible manner,
    - PATCH version when you make backward compatible bug fixes.

2. Make a tagged commit:
    ```shell
    git commit --message "Your commit message"
    git tag --annotate vMAJOR.MINOR.PATCH --message "vMAJOR.MINOR.PATCH"
    git push origin main --tags
    ```

3. Clean previous builds:
    ```shell
    rm -rf dist/*
    ```

4. Build the package:
    ```shell
    python3 -m build
    ```

5. Upload the package to TestPyPI and PyPI:
    ```shell
    # TestPyPI
    python3 -m twine upload --repository testpypi dist/*

    # PyPI
    python3 -m twine upload dist/*
    ```
