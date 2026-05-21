# Biomass Calculations

Small Python library for biomass calculations based on bundled species coefficients and lookup data.

The codebase is centered on one runtime module, `src/biomass.py`, plus one CSV dataset at `src/supplementary_data/Biomass_EquationParameters.csv`.


## What This Project Does

The library exposes three public functions:

- `getTreeBiomass(...)`
  - Calculates species-specific biomass for tree components such as wood, bark, branches, and foliage

- `getPhotoloadBiomass(...)`
  - Calculates biomass for supported Photoload plant codes

- `getDuffLitterBiomass(...)`
  - Returns duff/litter bulk density values or biomass estimates

All three functions support scalar inputs and vector-style sequence inputs.


## Codebase Shape

### Main files

- `src/biomass.py`
  - Main implementation
  - Loads the CSV data at import time
  - Defines constants, helper functions, and the public API

- `src/supplementary_data/Biomass_EquationParameters.csv`
  - Species metadata and coefficient dataset used by the runtime

- `tests/test_biomass.py`
  - `unittest` suite for the public API

- `docs/CODEBASE.md`
  - Longer architecture and codebase reference for this repo

### Architecture summary

At a high level:

1. Python imports `biomass`
2. `src/biomass.py` loads the bundled CSV into memory
3. callers invoke one of the public functions
4. inputs are validated and normalized
5. the calculation runs against in-memory coefficients and lookup tables


## Requirements

- Python 3.11+ recommended
- No third-party packages required for the current codebase

This project currently uses only the Python standard library.


## Installation

Because the code lives under `src/` and the repo does not currently include packaging metadata such as `pyproject.toml`, the simplest way to use it is to run with `src/` on `PYTHONPATH` or append `src/` to `sys.path`.

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

There are no external dependencies to install after that.

### Option 2: Use the repo directly

From the repository root:

Windows PowerShell:

```powershell
$env:PYTHONPATH = 'src'
```

macOS/Linux:

```bash
export PYTHONPATH=src
```


## Running the Code

### Quick import check

Windows PowerShell:

```powershell
$env:PYTHONPATH = 'src'
python -c "import biomass; print('loaded', len(biomass.BIOMASS_DATA))"
```

macOS/Linux:

```bash
PYTHONPATH=src python3 -c "import biomass; print('loaded', len(biomass.BIOMASS_DATA))"
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
import biomass

values = biomass.getTreeBiomass(
    spp=['PY', 'FDI'],
    decayclass=[1, 4],
    components=['wood', 'foliage'],
    dbh=[30.0, 40.0],
    height=[15.0, 20.0],
)

print(values)
```

### Example: Photoload biomass

```python
import biomass

value = biomass.getPhotoloadBiomass('AMAL', 50.0)
print(value)
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
import biomass

value = biomass.getDuffLitterBiomass(
    spp='PY',
    return_type='biomass',
    duff_depth=3.4,
    litter_depth=0.7,
)

print(value)
```


## Running Tests

From the repository root:

```powershell
python -m unittest discover -s tests -v
```

The tests currently adjust `sys.path` themselves so they can import `src/biomass.py`.


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

- a float for a single component scalar call
- a tuple for a multi-component scalar call
- a list for vectorized calls

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


## Data and Path Assumptions

Some behavior depends on the current repo layout:

- `src/biomass.py` loads its CSV relative to `__file__`
- the CSV must remain at `src/supplementary_data/Biomass_EquationParameters.csv`
- importing `biomass` performs file I/O immediately because the CSV is loaded at import time


## Notes

- Scalar and vector behavior share the same internal calculation path
- Invalid inputs generally raise `TypeError` or `ValueError`
- Invalid Photoload species emit a warning and return `0.0`


## Additional Reference

For a fuller architecture summary, see [docs/CODEBASE.md](docs/CODEBASE.md).
