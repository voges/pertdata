# Functionality for Handling Perturbation Datasets

Install the `pertdata` package:

```shell
pip3 install --requirement requirements.txt
pip3 install --editable .
```

Use the `PertDataset` class:

```python
from pertdata import PertDataset

pdata = PertDataset(
    name="replogle_2020_exp7", variant="preprocessed", dir_path="datasets"
)

print(pdata)
```

Please also see our [tutorial](tutorial.ipynb).
