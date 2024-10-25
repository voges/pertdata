# Contributing

## Package and Environment Management

We use [pip](https://pip.pypa.io) for package and environment management.

Here are some useful commands:

```sh
# Create a virtual environment.
python3 -m venv .venv

# Activate a virtual environment.
source .venv/bin/activate

# Install packages from a requirements.txt file.
pip3 install -r requirements.txt

# Install a package.
pip3 install <package>

# Update a requirements file.
pip3 freeze > requirements.txt

# Deactivate a virtual environment.
deactivate
```

## Linting & Formatting

We use [Ruff](https://github.com/astral-sh/ruff) to lint and format the code.
We recommend using the [Ruff extension for Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff).
