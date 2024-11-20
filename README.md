# Functionality for Handling Perturbation Datasets

Install the `pertdata` package:

```shell
pip3 install --requirement requirements.txt
pip3 install --editable .
```

Run the preprocessing:

```shell
pertdata-preprocess --help
```

Use the `PertDataset` class:

```python
from pertdata import PertDataset

ds = PertDataset(name="norman", variant="gears", dir_path="datasets")
print(ds)
```
