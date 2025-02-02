# Contributing

## Development

This repository provides a [development container](https://code.visualstudio.com/docs/devcontainers/containers) (or devcontainer).

In the devcontainer, the `pertdata` package is already installed (see [post_create.sh](.devcontainer/post_create.sh)).

The project is configured via [pyproject.toml](pyproject.toml).

## Distribution

1. Update the version in [`src/pertdata/version.py`](src/pertdata/version.py). Use [Semantic Versioning](https://semver.org).

2. Make a tagged commit:
    ```shell
    git add src/pertdata/version.py
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
    hatch build
    ```

5. Upload the package to TestPyPI and PyPI:
    ```shell
    # TestPyPI
    python3 -m twine upload --repository testpypi dist/*

    # PyPI
    python3 -m twine upload dist/*
    ```
