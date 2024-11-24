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

To build the package, execute the following command:
```python
python3 -m build
```

To upload the package to TestPyPI or PyPI, execute one of the following commands:
```python
# TestPyPI
python3 -m twine upload --repository testpypi dist/*

# PyPI
python3 -m twine upload dist/*
```

To install the package from TestPyPI or PyPI, execute one of the following commands:
```python
# TestPyPI
pip3 install --index-url https://test.pypi.org/simple/ --no-deps pertdata

# PyPI
pip3 install pertdata
```
