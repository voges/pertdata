# Functionality for Handling Perturbation Datasets

Install the `pertdata` package:

```sh
bash install.sh
```

Run the preprocessing:

```sh
pertdata-preprocess --help
```

Use the `PertDataset` class:

```python
from pertdata import PertDataset

ds = PertDataset(name="norman", variant="gears", dir_path="datasets")
print(ds)
```
