__author__ = ['Gregory A. Greene, map.n.trowel@gmail.com']

import csv
import math
import numpy as np
import warnings
from collections.abc import Sequence
from numbers import Real
from pathlib import Path
from typing import Optional, Union


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


def _as_scalar_or_vector(results: list, is_vectorized: bool):
    return results if is_vectorized else results[0]


def _any_array(*args) -> bool:
    """Return True if any argument is a np.ndarray."""
    return any(isinstance(a, np.ndarray) for a in args)


def _to_1d_array(value) -> np.ndarray:
    """Coerce a scalar or array to a 1-D ndarray."""
    return np.atleast_1d(np.asarray(value))


def _broadcast_arguments(**kwargs):
    lengths = {len(value) for value in kwargs.values() if _is_sequence(value)}
    if not lengths:
        return False, {name: [value] for name, value in kwargs.items()}
    if len(lengths) != 1:
        raise ValueError('Vector inputs must all be the same length')

    vector_length = lengths.pop()
    broadcast = {}
    for name, value in kwargs.items():
        if _is_sequence(value):
            broadcast[name] = list(value)
        else:
            broadcast[name] = [value] * vector_length
    return True, broadcast


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


def _get_decay_vector(species: str, decayclass: int) -> tuple[float, ...]:
    decay_lookup = SOFTWOOD_DECAY if species in SOFTWOOD_SPECIES else HARDWOOD_DECAY
    try:
        return decay_lookup[decayclass]
    except KeyError as exc:
        tree_type = 'softwood' if species in SOFTWOOD_SPECIES else 'hardwood'
        raise ValueError(f'Invalid decay class "{decayclass}" for {tree_type} species "{species}"') from exc


def _get_species_row(species: str) -> dict[str, Union[str, float]]:
    try:
        return BIOMASS_DATA[species]
    except KeyError as exc:
        raise ValueError(f'Unknown species code: "{species}"') from exc


def _is_finite_number(value: object) -> bool:
    return isinstance(value, Real) and math.isfinite(value)


def _is_sequence(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _load_biomass_data(path: Path) -> dict[str, dict[str, Union[str, float]]]:
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        rows = {}
        for row in reader:
            species = row['Species']
            parsed = {}
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
    elif _is_sequence(components):
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
    if _is_sequence(depth):
        normalized = []
        for value in depth:
            if not isinstance(value, Real):
                raise TypeError(f'"{name}" values must be numeric and expressed in cm')
            normalized.append(value / 100)
        return normalized
    if not isinstance(depth, Real):
        raise TypeError(f'"{name}" must be numeric and expressed in cm')
    return depth / 100.0


def _resolve_photoload_height(pl_code: str, height: Optional[float]) -> Optional[float]:
    if height == 0:
        return 0.0
    if height is None or not _is_finite_number(height):
        return PHOTOLOAD_DEFAULT_HEIGHTS.get(pl_code)
    return height


def _validate_species_mix(
    spp: Union[str, Sequence[str]],
    pct_list: Optional[Sequence[float]],
) -> tuple[list[str], Optional[list[float]]]:
    if isinstance(spp, str):
        return [spp], None
    if not _is_sequence(spp):
        raise TypeError('"spp" must be a string or a sequence of species codes')
    if not all(isinstance(species, str) for species in spp):
        raise TypeError('Each species code in "spp" must be a string')
    if pct_list is None:
        raise ValueError('"pct_list" is required when "spp" contains multiple species')
    if not _is_sequence(pct_list):
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
    duff_depth: Optional[Union[int, float, Sequence[Union[int, float]]]] = None,
    litter_depth: Optional[Union[int, float, Sequence[Union[int, float]]]] = None,
) -> Union[float, tuple[float, float], list[float], list[tuple[float, float]]]:
    """
    Return duff/litter bulk density values or calculate duff/litter biomass.

    Species mixes are provided through ``spp`` + ``pct_list``. ``duff_depth`` and ``litter_depth`` may be
    scalars or same-length vectors expressed in cm.
    """
    species_list, pct_values = _validate_species_mix(spp, pct_list)
    if return_type not in {'bulk_density', 'biomass'}:
        raise ValueError('"return_type" must be either "bulk_density" or "biomass"')

    duff_bd, litter_bd = _get_bulk_density(species_list, pct_values)
    if return_type == 'bulk_density':
        return duff_bd, litter_bd

    normalized_duff = _normalize_depth(duff_depth, 'duff_depth')
    normalized_litter = _normalize_depth(litter_depth, 'litter_depth')
    if normalized_duff is None and normalized_litter is None:
        raise ValueError('At least one of "duff_depth" or "litter_depth" must be provided for biomass')

    is_vectorized, args = _broadcast_arguments(duff_depth=normalized_duff, litter_depth=normalized_litter)
    results = []
    for duff_value, litter_value in zip(args['duff_depth'], args['litter_depth']):
        if duff_value is not None and litter_value is not None:
            results.append((duff_bd * duff_value, litter_bd * litter_value))
        elif duff_value is not None:
            results.append(duff_bd * duff_value)
        else:
            results.append(litter_bd * litter_value)

    return _as_scalar_or_vector(results, is_vectorized)


def getPhotoloadBiomass(
    pl_code: Union[str, Sequence[str]],
    pct_cvr: Union[float, Sequence[float]],
    height: Optional[Union[float, Sequence[Optional[float]]]] = None,
) -> Union[float, list[float]]:
    """
    Estimate Photoload biomass in kg/m2.

    Scalar and vector inputs are supported for ``pl_code``, ``pct_cvr``, and ``height``.
    Vector inputs must all have the same length.
    """
    is_vectorized, args = _broadcast_arguments(pl_code=pl_code, pct_cvr=pct_cvr, height=height)

    results = []
    for species, cover, species_height in zip(args['pl_code'], args['pct_cvr'], args['height']):
        if cover == 0 or species_height == 0:
            results.append(0.0)
            continue

        resolved_height = _resolve_photoload_height(species, species_height)
        if resolved_height is None:
            warnings.warn(f'Photoload species code "{species}" is invalid. Biomass returned as 0.')
            results.append(0.0)
            continue

        results.append(_calculate_photoload_biomass(species, cover, resolved_height))

    return _as_scalar_or_vector(results, is_vectorized)


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
    return_array = _any_array(spp, decayclass, dbh, height)

    spp_arr = _to_1d_array(spp)
    dc_arr = _to_1d_array(decayclass).astype(int)
    dbh_arr = _to_1d_array(dbh).astype(float)
    ht_arr = _to_1d_array(height).astype(float) if height is not None else None

    # Validate that explicitly-provided ndarray inputs all share the same length.
    # Scalar inputs (non-ndarray) are allowed to broadcast freely.
    explicit_array_inputs = [x for x in (spp, decayclass, dbh, height) if isinstance(x, np.ndarray)]
    explicit_lengths = {arr.shape[0] for arr in explicit_array_inputs}
    if len(explicit_lengths) > 1:
        raise ValueError('Vector inputs must all be the same length')

    # Determine n from multi-element arrays; broadcast length-1 arrays (from scalars)
    candidate_arrs = [spp_arr, dc_arr, dbh_arr]
    if ht_arr is not None:
        candidate_arrs.append(ht_arr)
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
