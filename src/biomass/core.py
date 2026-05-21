__author__ = ['Gregory A. Greene, map.n.trowel@gmail.com']

import csv
import math
import numpy as np
import numpy.typing as npt
import warnings
from numbers import Real
from pathlib import Path
from typing import Sequence


DATA_PATH = Path(__file__).resolve().parent / 'supplementary_data' / 'Biomass_EquationParameters.csv'
TREE_COMPONENTS = ('wood', 'bark', 'branches', 'foliage')
TREE_COMPONENT_INDEX = {component: index for index, component in enumerate(TREE_COMPONENTS)}
SOFTWOOD_DECAY = {
    1: (1, 1, 1, 1),
    2: (1, 1, 1, 1),
    3: (1, 1, 1, 0),
    4: (0.95, 0.75, 0.5, 0),
    5: (0.80, 0.5, 0, 0),
    6: (0.66, 0.25, 0, 0),
    7: (0.5, 0.1, 0, 0),
    8: (0.33, 0.05, 0, 0),
    9: (0.05, 0, 0, 0),
}
HARDWOOD_DECAY = {
    1: (1, 1, 1, 1),
    2: (1, 1, 1, 1),
    3: (1, 1, 1, 0),
    4: (0.95, 0.75, 0.5, 0),
    5: (0.66, 0.25, 0, 0),
    6: (0.05, 0, 0, 0),
}
PHOTOLOAD_DEFAULT_HEIGHTS = {
    'AMAL': 35.56,
    'BERE': 10.16,
    'MARE': 10.16,
    'SYAL': 35.56,
    'VAGL': 35.56,
    'VASC': 17.78,
    'ARLA': 30.48,
    'CARU': 20.3,
    'FESC': 30.5,
    'XETE': 25.4,
}


def _any_array(*args) -> bool:
    """Return True if any argument is a np.ndarray."""
    return any(isinstance(a, np.ndarray) for a in args)


def _to_1d_array(value) -> npt.NDArray:
    """Coerce a scalar or array to a 1-D ndarray."""
    return np.atleast_1d(np.asarray(value))


def _calculate_photoload_biomass(pl_code: str, pct_cvr, height):
    if pl_code == 'AMAL':
        return (height / 35.56) * 0.0148 * np.exp(0.0454 * pct_cvr)
    if pl_code in {'BERE', 'MARE'}:
        return (height / 10.16) * 0.0013 * pct_cvr
    if pl_code == 'SYAL':
        return (height / 35.56) * 0.0107 * np.exp(0.0442 * pct_cvr)
    if pl_code == 'VAGL':
        return (height / 35.56) * 0.0052 * pct_cvr
    if pl_code == 'VASC':
        return (height / 17.78) * 0.0135 * np.exp(0.035 * pct_cvr)
    if pl_code == 'ARLA':
        return (height / 30.48) * ((0.000008 * pct_cvr ** 2) + (0.0006 * pct_cvr) + 0.0022)
    if pl_code == 'CARU':
        return (height / 20.3) * 0.0194 * np.exp(0.0282 * pct_cvr)
    if pl_code == 'FESC':
        return (height / 30.5) * 0.0188 * np.exp(0.0298 * pct_cvr)
    if pl_code == 'XETE':
        return (height / 25.4) * 0.0154 * np.exp(0.0444 * pct_cvr)
    warnings.warn(f'Photoload species code "{pl_code}" is invalid. Biomass returned as 0.')
    return 0.0


def _get_bulk_density(spp: list[str], pct_list: Optional[list[float]]) -> tuple[float, float]:
    if pct_list is None:
        row = _get_species_row(spp[0])
        return float(row['DUFF_BD']), float(row['LITTER_BD'])

    duff_bd = 0.0
    litter_bd = 0.0
    for species, pct in zip(spp, pct_list):
        row = _get_species_row(species)
        duff_bd += row['DUFF_BD'] * pct / 100
        litter_bd += row['LITTER_BD'] * pct / 100
    return duff_bd, litter_bd


def _get_species_row(species: str) -> dict[str, Union[str, float]]:
    try:
        return BIOMASS_DATA[species]
    except KeyError as exc:
        raise ValueError(f'Unknown species code: "{species}"') from exc


def _load_biomass_data(path: Path) -> dict[str, dict[str, Union[str, float]]]:
    rows: dict[str, dict[str, str | float]] = {}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            species = row['Species']
            parsed: dict[str, str | float] = {}
            for key, value in row.items():
                if key in {
                    'Species',
                    'Species_Label',
                    'FOFEM_sppCD',
                    'GCBM_GenSppCD',
                    'TreeType',
                    'DUFF_LITTER_BIOMASS_Ref',
                    'TREE_BIOMASS_Ref',
                }:
                    parsed[key] = value
                else:
                    parsed[key] = float(value)
            rows[species] = parsed
    return rows


def _normalize_components(components: Union[str, Sequence[str]]) -> list[str]:
    if isinstance(components, str):
        component_list = [components]
    elif isinstance(components, Sequence) and not isinstance(components, (str, bytes, bytearray)):
        component_list = list(components)
    else:
        raise TypeError('"components" must be a string or a sequence of strings')

    invalid_components = [component for component in component_list if component not in TREE_COMPONENTS]
    if invalid_components:
        raise ValueError(f'Invalid tree components: {invalid_components}')
    return component_list


def _normalize_depth(depth: Optional[Union[int, float, np.ndarray]], name: str):
    if depth is None:
        return None
    if isinstance(depth, np.ndarray):
        if not np.issubdtype(depth.dtype, np.number):
            raise TypeError(f'"{name}" values must be numeric and expressed in cm')
        return depth / 100.0
    if not isinstance(depth, Real):
        raise TypeError(f'"{name}" must be numeric and expressed in cm')
    return depth / 100.0


def _validate_species_mix(
    spp: Union[str, Sequence[str]],
    pct_list: Optional[Sequence[float]],
) -> tuple[list[str], Optional[list[float]]]:
    if isinstance(spp, str):
        return [spp], None
    if not (isinstance(spp, Sequence) and not isinstance(spp, (str, bytes, bytearray))):
        raise TypeError('"spp" must be a string or a sequence of species codes')
    if not all(isinstance(species, str) for species in spp):
        raise TypeError('Each species code in "spp" must be a string')
    if pct_list is None:
        raise ValueError('"pct_list" is required when "spp" contains multiple species')
    if not (isinstance(pct_list, Sequence) and not isinstance(pct_list, (str, bytes, bytearray))):
        raise TypeError('"pct_list" must be a sequence of numeric percentages')
    if len(spp) != len(pct_list):
        raise ValueError('"spp" and "pct_list" must be the same length')
    if not all(isinstance(pct, Real) for pct in pct_list):
        raise TypeError('Each value in "pct_list" must be numeric')
    if not math.isclose(sum(pct_list), 100.0, abs_tol=0.5):
        raise ValueError('The values in "pct_list" must sum to 100')
    return list(spp), list(pct_list)


BIOMASS_DATA = _load_biomass_data(DATA_PATH)
SOFTWOOD_SPECIES = {species for species, row in BIOMASS_DATA.items() if row['TreeType'] == 'softwood'}


def getDuffLitterBiomass(
    spp: Union[str, Sequence[str]],
    pct_list: Optional[Sequence[float]] = None,
    return_type: str = 'bulk_density',
    duff_depth: Optional[Union[int, float, np.ndarray]] = None,
    litter_depth: Optional[Union[int, float, np.ndarray]] = None,
) -> Union[float, tuple[float, float], np.ndarray, tuple[np.ndarray, np.ndarray]]:
    """
    Return duff/litter bulk density values or calculate duff/litter biomass.

    Parameters
    ----------
    spp : str or sequence of str
        Species code or sequence of species codes defining the stand composition.
        When a sequence is given, ``pct_list`` is required.
    pct_list : sequence of float, optional
        Percentage of each species in ``spp``. Values must sum to 100 (±0.5).
        Required when ``spp`` contains more than one species.
    return_type : {'bulk_density', 'biomass'}
        'bulk_density' returns the weighted duff and litter bulk density values
        in kg/m³. 'biomass' multiplies bulk density by the supplied depths to
        return biomass in kg/m².
    duff_depth : int, float, or np.ndarray, optional
        Duff layer depth in cm. Scalar or 1-D array. Required when
        ``return_type='biomass'`` and ``litter_depth`` is not given.
    litter_depth : int, float, or np.ndarray, optional
        Litter layer depth in cm. Scalar or 1-D array. Required when
        ``return_type='biomass'`` and ``duff_depth`` is not given.

    Returns
    -------
    tuple of float
        ``(duff_bulk_density, litter_bulk_density)`` in kg/m³ when
        ``return_type='bulk_density'``.
    float
        Biomass in kg/m² when scalar depth inputs and only one depth is given.
    tuple of float
        ``(duff_biomass, litter_biomass)`` in kg/m² when scalar depth inputs
        and both depths are given.
    np.ndarray
        Biomass array in kg/m² when any depth input is np.ndarray and only
        one depth is given.
    tuple of np.ndarray
        ``(duff_biomass, litter_biomass)`` arrays in kg/m² when any depth input
        is np.ndarray and both depths are given.

    Raises
    ------
    ValueError
        If ``return_type`` is invalid, if ``pct_list`` is missing for a
        multi-species ``spp``, if ``pct_list`` does not sum to 100, or if no
        depth is provided when ``return_type='biomass'``.
    TypeError
        If ``spp`` or ``pct_list`` elements are the wrong type, or if depth
        values are non-numeric.
    """
    species_list, pct_values = _validate_species_mix(spp, pct_list)
    if return_type not in {'bulk_density', 'biomass'}:
        raise ValueError('"return_type" must be either "bulk_density" or "biomass"')

    duff_bd, litter_bd = _get_bulk_density(species_list, pct_values)
    if return_type == 'bulk_density':
        return duff_bd, litter_bd

    return_array = _any_array(duff_depth, litter_depth)

    normalized_duff = _normalize_depth(duff_depth, 'duff_depth')
    normalized_litter = _normalize_depth(litter_depth, 'litter_depth')

    if isinstance(normalized_duff, np.ndarray) and isinstance(normalized_litter, np.ndarray):
        if normalized_duff.shape != normalized_litter.shape:
            raise ValueError('"duff_depth" and "litter_depth" must have the same length')

    if normalized_duff is None and normalized_litter is None:
        raise ValueError('At least one of "duff_depth" or "litter_depth" must be provided for biomass')

    if normalized_duff is not None and normalized_litter is not None:
        duff_arr = np.atleast_1d(np.asarray(duff_bd * normalized_duff, dtype=float))
        litter_arr = np.atleast_1d(np.asarray(litter_bd * normalized_litter, dtype=float))
        if return_array:
            n = max(duff_arr.shape[0], litter_arr.shape[0])
            if duff_arr.shape[0] == 1 and n > 1:
                duff_arr = np.repeat(duff_arr, n)
            if litter_arr.shape[0] == 1 and n > 1:
                litter_arr = np.repeat(litter_arr, n)
            return duff_arr, litter_arr
        return float(duff_arr[0]), float(litter_arr[0])
    elif normalized_duff is not None:
        result = np.asarray(duff_bd * normalized_duff)
        return np.atleast_1d(result) if return_array else float(result)
    else:
        result = np.asarray(litter_bd * normalized_litter)
        return np.atleast_1d(result) if return_array else float(result)


def getPhotoloadBiomass(
    pl_code: Union[str, np.ndarray],
    pct_cvr: Union[float, np.ndarray],
    height: Optional[Union[float, np.ndarray]] = None,
) -> Union[float, np.ndarray]:
    """
    Estimate Photoload biomass in kg/m².

    Parameters
    ----------
    pl_code : str or np.ndarray
        Photoload species code or 1-D array of species codes.
        Supported codes: 'AMAL', 'BERE', 'MARE', 'SYAL', 'VAGL', 'VASC',
        'ARLA', 'CARU', 'FESC', 'XETE'. Invalid codes produce a warning and
        return 0.0 for the affected elements.
    pct_cvr : float or np.ndarray
        Percent cover (0-100). Scalar or 1-D array.
    height : float or np.ndarray, optional
        Canopy height in cm. Scalar or 1-D array. Pass ``None`` to use the
        default height for all elements. Within an array, use ``np.nan`` to
        use the species default height for individual elements.

    Returns
    -------
    float
        When all inputs are scalars.
    np.ndarray
        When any input is np.ndarray.

    Raises
    ------
    ValueError
        If array inputs have mismatched lengths or if an ndarray is not 1-D.
    TypeError
        If plain sequences (list, tuple) are passed instead of np.ndarray.
    """
    # Validate: reject plain sequences and 0-D arrays
    for _name, _val in [('pl_code', pl_code), ('pct_cvr', pct_cvr)] + (
        [('height', height)] if height is not None else []
    ):
        if isinstance(_val, Sequence) and not isinstance(_val, (str, bytes, bytearray)):
            raise TypeError(
                f'"{_name}" does not accept plain sequences; use np.array() to convert'
            )
        if isinstance(_val, np.ndarray) and _val.ndim != 1:
            raise ValueError(f'"{_name}" must be a 1-D array, got shape {_val.shape}')

    return_array = _any_array(pl_code, pct_cvr, height)

    pl_code_arr = _to_1d_array(pl_code)
    pct_cvr_arr = _to_1d_array(pct_cvr).astype(float)

    # height=None → all defaults (represented internally as nan)
    if height is None:
        ht_arr = np.full(pl_code_arr.shape[0], np.nan)
    else:
        ht_arr = _to_1d_array(height).astype(float)

    # Raise early if two caller-supplied ndarrays have different lengths
    explicit_ndarray_lengths = {
        arr.shape[0]
        for val, arr in [(pl_code, pl_code_arr), (pct_cvr, pct_cvr_arr)]
        + ([(height, ht_arr)] if isinstance(height, np.ndarray) else [])
        if isinstance(val, np.ndarray)
    }
    if len(explicit_ndarray_lengths) > 1:
        raise ValueError('Vector inputs must all be the same length')

    # Broadcast length-1 arrays to match multi-element arrays
    multi_lengths = {arr.shape[0] for arr in [pl_code_arr, pct_cvr_arr, ht_arr] if arr.shape[0] > 1}
    if len(multi_lengths) > 1:
        raise ValueError('Vector inputs must all be the same length')
    n = multi_lengths.pop() if multi_lengths else 1

    pl_code_arr = np.repeat(pl_code_arr, n) if pl_code_arr.shape[0] == 1 else pl_code_arr
    pct_cvr_arr = np.repeat(pct_cvr_arr, n) if pct_cvr_arr.shape[0] == 1 else pct_cvr_arr
    ht_arr = np.repeat(ht_arr, n) if ht_arr.shape[0] == 1 else ht_arr

    output = np.zeros(n)

    for code in np.unique(pl_code_arr):
        code_str = str(code)
        code_mask = (pl_code_arr == code)
        default_height = PHOTOLOAD_DEFAULT_HEIGHTS.get(code_str)

        h_slice = ht_arr[code_mask]
        pct_slice = pct_cvr_arr[code_mask]

        # Replace nan with the code's default height; if no default, nan remains
        if default_height is not None:
            resolved = np.where(np.isnan(h_slice), default_height, h_slice)
        else:
            resolved = h_slice.copy()

        if np.any(np.isnan(resolved)):
            warnings.warn(f'Photoload species code "{code_str}" is invalid. Biomass returned as 0.')
            resolved = np.where(np.isnan(resolved), 0.0, resolved)

        zero_mask = (pct_slice == 0) | (resolved == 0)
        active = ~zero_mask

        if np.any(active):
            result_slice = np.zeros(int(np.sum(code_mask)))
            result_slice[active] = _calculate_photoload_biomass(
                code_str, pct_slice[active], resolved[active]
            )
            output[code_mask] = result_slice

    return output if return_array else float(output[0])


def getTreeBiomass(
    spp: Union[str, np.ndarray],
    decayclass: Union[int, np.ndarray],
    components: Union[str, Sequence[str]],
    dbh: Union[float, np.ndarray],
    height: Optional[Union[float, np.ndarray]] = None,
) -> Union[float, tuple[float, ...], np.ndarray, tuple[np.ndarray, ...]]:
    """
    Return species-specific biomass values for tree components using Canadian National Biomass equations.

    Parameters
    ----------
    spp : str or np.ndarray
        Species code or 1-D array of species codes.
    decayclass : int or np.ndarray
        Decay class integer or 1-D array of decay class integers.
        Valid range is 1-9 for softwood species and 1-6 for hardwood species.
    components : str or sequence of str
        One or more of 'wood', 'bark', 'branches', 'foliage'.
    dbh : float or np.ndarray
        Diameter at breast height in cm. Scalar or 1-D array.
    height : float or np.ndarray, optional
        Tree height in m. Scalar or 1-D array. When omitted, the DBH-only
        allometric equation is used for all trees in the call.

    Returns
    -------
    float
        When all inputs are scalars and a single component is requested.
    tuple of float
        When all inputs are scalars and multiple components are requested.
    np.ndarray
        When any input is np.ndarray and a single component is requested.
    tuple of np.ndarray
        When any input is np.ndarray and multiple components are requested.

    Raises
    ------
    ValueError
        If an unknown species code is given, if the decay class is out of range
        for the species type, or if array inputs have mismatched lengths.
    TypeError
        If components is not a string or sequence of strings.
    """
    component_list = _normalize_components(components)

    # Reject non-1-D ndarrays (e.g. 0-D scalars wrapped in np.array) before any conversion.
    for _name, _val in [('spp', spp), ('decayclass', decayclass), ('dbh', dbh)] + (
        [('height', height)] if isinstance(height, np.ndarray) else []
    ):
        if isinstance(_val, np.ndarray) and _val.ndim != 1:
            raise ValueError(f'"{_name}" must be a 1-D array, got shape {_val.shape}')

    # Reject plain Python sequences — callers must use np.array() for vector inputs.
    for _name, _val in [('spp', spp), ('decayclass', decayclass), ('dbh', dbh)] + (
        [('height', height)] if height is not None else []
    ):
        if isinstance(_val, Sequence) and not isinstance(_val, (str, bytes, bytearray)):
            raise TypeError(
                f'"{_name}" does not accept plain sequences; use np.array() to convert'
            )

    return_array = _any_array(spp, decayclass, dbh, height)

    spp_arr = _to_1d_array(spp)
    dc_arr = _to_1d_array(decayclass).astype(int)
    dbh_arr = _to_1d_array(dbh).astype(float)
    ht_arr = _to_1d_array(height).astype(float) if height is not None else None

    # Validate that explicitly-provided ndarray inputs all share the same length.
    # Scalar inputs (non-ndarray) are allowed to broadcast freely.
    explicit_pairs = [('spp', spp, spp_arr), ('decayclass', decayclass, dc_arr), ('dbh', dbh, dbh_arr)]
    if ht_arr is not None:
        explicit_pairs.append(('height', height, ht_arr))
    explicit_lengths = {arr.shape[0] for _, orig, arr in explicit_pairs if isinstance(orig, np.ndarray)}
    if len(explicit_lengths) > 1:
        raise ValueError('Vector inputs must all be the same length')

    # Determine n from multi-element arrays; broadcast length-1 arrays (from scalars).
    candidate_arrs = [spp_arr, dc_arr, dbh_arr] + ([ht_arr] if ht_arr is not None else [])
    multi_lengths = {arr.shape[0] for arr in candidate_arrs if arr.shape[0] > 1}
    n = multi_lengths.pop() if multi_lengths else 1

    spp_arr = np.repeat(spp_arr, n) if spp_arr.shape[0] == 1 else spp_arr
    dc_arr = np.repeat(dc_arr, n) if dc_arr.shape[0] == 1 else dc_arr
    dbh_arr = np.repeat(dbh_arr, n) if dbh_arr.shape[0] == 1 else dbh_arr
    if ht_arr is not None:
        ht_arr = np.repeat(ht_arr, n) if ht_arr.shape[0] == 1 else ht_arr

    outputs = {comp: np.zeros(n) for comp in component_list}

    for species in np.unique(spp_arr):
        species_str = str(species)
        species_mask = (spp_arr == species)
        row = _get_species_row(species_str)
        decay_lookup = SOFTWOOD_DECAY if species_str in SOFTWOOD_SPECIES else HARDWOOD_DECAY

        for dc in np.unique(dc_arr[species_mask]):
            dc_mask = species_mask & (dc_arr == dc)
            try:
                decay_vector = decay_lookup[int(dc)]
            except KeyError as exc:
                tree_type = 'softwood' if species_str in SOFTWOOD_SPECIES else 'hardwood'
                raise ValueError(
                    f'Invalid decay class "{dc}" for {tree_type} species "{species_str}"'
                ) from exc

            for component in component_list:
                decay_factor = decay_vector[TREE_COMPONENT_INDEX[component]]
                if ht_arr is None:
                    outputs[component][dc_mask] = (
                        decay_factor
                        * row[f'BIOMASS_DBH_B{component}1']
                        * np.power(dbh_arr[dc_mask], row[f'BIOMASS_DBH_B{component}2'])
                    )
                else:
                    outputs[component][dc_mask] = (
                        decay_factor
                        * row[f'BIOMASS_DBH-HT_B{component}1']
                        * np.power(dbh_arr[dc_mask], row[f'BIOMASS_DBH-HT_B{component}2'])
                        * np.power(ht_arr[dc_mask], row[f'BIOMASS_DBH-HT_B{component}3'])
                    )

    if isinstance(components, str):
        result = outputs[component_list[0]]
        return result if return_array else float(result[0])
    else:
        arrays = tuple(outputs[comp] for comp in component_list)
        return arrays if return_array else tuple(float(arr[0]) for arr in arrays)
