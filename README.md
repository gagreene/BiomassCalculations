# Biomass Calculations

Small Python library for biomass calculations based on bundled species coefficients and lookup data.

The codebase is centered on one runtime package, `src/biomass/`, with implementation in `src/biomass/core.py` and one CSV dataset at `src/biomass/supplementary_data/Biomass_EquationParameters.csv`.


## What This Project Does

The library exposes three public functions:

- `getTreeBiomass(...)`
  - Calculates species-specific biomass for tree components such as wood, bark, branches, and foliage
  - Uses the Canadian National Biomass Equations from Lambert et al. (2005) and Ung et al. (2008)

- `getPhotoloadBiomass(...)`
  - Calculates biomass for supported Photoload plant codes

- `getDuffLitterBiomass(...)`
  - Returns duff/litter bulk density values or biomass estimates

All three functions support scalar inputs and `np.ndarray` inputs.


## Scientific Basis

The tree biomass calculations in this package are based on the Canadian National Biomass Equations, using equations from:

- Lambert, M.-C., Ung, C.-H., and Raulier, F. (2005). *Canadian national tree aboveground biomass equations.*
- Ung, C.-H., Bernier, P., Guo, X.-J., Lambert, M.-C., and Regniere, J. (2008). *Canadian national biomass equations.*

Those references provide the basis for the tree component biomass coefficients bundled in the package dataset.


## Codebase Shape

### Main files

- `src/biomass/__init__.py`
  - Package entry point
  - Re-exports the public API from `src/biomass/core.py`

- `src/biomass/core.py`
  - Main implementation
  - Loads the CSV data at import time
  - Defines constants, helper functions, and the public API logic

- `src/biomass/supplementary_data/Biomass_EquationParameters.csv`
  - Species metadata and coefficient dataset used by the runtime

- `tests/test_biomass.py`
  - `unittest` suite for the public API

- `docs/CODEBASE.md`
  - Longer architecture and codebase reference for this repo

### Architecture summary

At a high level:

1. Python imports `biomass`
2. `src/biomass/core.py` loads the bundled CSV into memory
3. callers invoke one of the public functions
4. inputs are validated and normalized
5. the calculation runs against in-memory coefficients and lookup tables


## Requirements

- Python 3.11+
- [NumPy](https://numpy.org/) >= 1.24

This project uses only the Python standard library and NumPy. No other third-party packages are required.


## Installation

This repo now includes `pyproject.toml`, so it can be installed as a standard Python package.

### Option 1: Create and use a virtual environment

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the package in editable mode from the repository root:

```powershell
python -m pip install -e .
```

For a regular local install:

```powershell
python -m pip install .
```

### Option 2: Use the repo directly without installation

If you are working directly from the repo, you can still run with `src/` on `PYTHONPATH`.


## Running the Code

### Quick import check

Windows PowerShell:

```powershell
python -c "import biomass; print('loaded', biomass.__version__, len(biomass.BIOMASS_DATA))"
```

macOS/Linux:

```bash
PYTHONPATH=src python3 -c "import biomass; print('loaded', biomass.__version__, len(biomass.BIOMASS_DATA))"
```

### Example: tree biomass

```python
import biomass

value = biomass.getTreeBiomass(
    spp='PY',
    decayclass=1,
    components='wood',
    dbh=30.0,
)

print(value)
```

### Example: tree biomass with multiple records

```python
import numpy as np
import biomass

values = biomass.getTreeBiomass(
    spp=np.array(['PY', 'FDI']),
    decayclass=np.array([1, 4]),
    components=['wood', 'foliage'],
    dbh=np.array([30.0, 40.0]),
    height=np.array([15.0, 20.0]),
)

print(values)  # tuple of two np.ndarray (wood values, foliage values)
```

### Example: Photoload biomass

```python
import numpy as np
import biomass

# Scalar call
value = biomass.getPhotoloadBiomass('AMAL', 50.0)
print(value)  # float

# Array call — pass np.nan in the height array to use the species default height
values = biomass.getPhotoloadBiomass(
    np.array(['AMAL', 'VAGL']),
    np.array([50.0, 30.0]),
    np.array([35.56, np.nan]),  # np.nan → uses VAGL default height
)
print(values)  # np.ndarray
```

### Example: duff/litter bulk density

```python
import biomass

value = biomass.getDuffLitterBiomass(
    spp=['PY', 'FDI', 'PLI'],
    pct_list=[60.0, 30.0, 10.0],
    return_type='bulk_density',
)

print(value)
```

### Example: duff/litter biomass

```python
import numpy as np
import biomass

# Scalar call
value = biomass.getDuffLitterBiomass(
    spp='PY',
    return_type='biomass',
    duff_depth=3.4,
    litter_depth=0.7,
)
print(value)  # (duff_biomass_float, litter_biomass_float)

# Array call — pass depth arrays to compute biomass for multiple plots
values = biomass.getDuffLitterBiomass(
    spp='PY',
    return_type='biomass',
    duff_depth=np.array([3.4, 5.0, 2.1]),
    litter_depth=np.array([0.7, 1.2, 0.9]),
)
print(values)  # (duff_ndarray, litter_ndarray)
```


## Running Tests

From the repository root:

```powershell
python -m unittest discover -s tests -v
```

The tests currently add `src/` to `sys.path` so they can run directly from the repo without requiring installation first.


## Packaging

### PyPI

This repo is prepared for Python packaging with:

- `pyproject.toml`
- a standard `src/` package layout
- bundled package data configuration for the CSV file

Build artifacts locally with:

```powershell
python -m build
```

### Conda

This repo also includes a starter Conda recipe:

- `conda-recipe/meta.yaml`

Build it with your Conda build tooling of choice, for example:

```powershell
conda build conda-recipe
```


## Public API

### `getTreeBiomass`

Calculates tree biomass for one or more species records.

Inputs:

- `spp`
- `decayclass`
- `components`
- `dbh`
- optional `height`

Returns:

- a `float` for a single component scalar call
- a `tuple[float, ...]` for a multi-component scalar call
- a `np.ndarray` for a single component array call
- a `tuple[np.ndarray, ...]` for a multi-component array call

### `getPhotoloadBiomass`

Calculates Photoload biomass for one or more records.

Inputs:

- `pl_code`
- `pct_cvr`
- optional `height`

### `getDuffLitterBiomass`

Returns bulk density or biomass for duff/litter.

Inputs:

- `spp`
- optional `pct_list`
- `return_type`
- optional `duff_depth`
- optional `litter_depth`


## Array vs Scalar Behavior

All three public functions support both scalar and array inputs:

- Pass plain Python scalars (`str`, `int`, `float`) to get scalar return values.
- Pass `np.ndarray` for any argument to get `np.ndarray` return values.
- Scalar arguments broadcast automatically when other arguments are arrays.

```python
import numpy as np
import biomass

# Scalar spp and decayclass broadcast across array dbh
result = biomass.getTreeBiomass('PY', 1, 'wood', np.array([20.0, 30.0, 40.0]))
# result is a np.ndarray of shape (3,)
```

Plain Python lists are **not** supported as vector inputs. Use `np.array(...)` to convert lists before calling.


## Data and Path Assumptions

Some behavior depends on the current repo layout:

- `src/biomass/core.py` loads its CSV relative to `__file__`
- the CSV must remain at `src/biomass/supplementary_data/Biomass_EquationParameters.csv`
- importing `biomass` performs file I/O immediately because the CSV is loaded at import time


## Notes

- Scalar and array inputs share the same internal calculation path via NumPy masked operations
- Invalid inputs generally raise `TypeError` or `ValueError`
- Invalid Photoload species emit a warning and return `0.0`
- Use `np.nan` within a height array to apply the species default height at individual positions


## Additional Reference

For a fuller architecture summary, see [docs/CODEBASE.md](docs/CODEBASE.md).
