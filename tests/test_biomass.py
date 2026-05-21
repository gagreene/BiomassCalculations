import math
import numpy as np
from pathlib import Path
import sys
import unittest
import warnings

PROJ_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ_ROOT / 'src'))

import biomass


class GetTreeBiomassTests(unittest.TestCase):
    """Behavioral tests for ``biomass.get_tree_biomass``."""

    def test_scalar_dbh_only_uses_dbh_coefficients(self):
        """Use DBH-only coefficients when ``height`` is omitted."""
        result = biomass.get_tree_biomass('PY', 1, 'wood', 30.0)
        expected = 0.0564 * math.pow(30.0, 2.4465)
        self.assertAlmostEqual(result, expected, places=9)

    def test_scalar_input_returns_float(self):
        """Return a float for scalar input and single component."""
        result = biomass.get_tree_biomass('PY', 1, 'wood', 30.0)
        self.assertIsInstance(result, float)
        expected = 0.0564 * math.pow(30.0, 2.4465)
        self.assertAlmostEqual(result, expected, places=9)

    def test_array_single_component_returns_ndarray(self):
        """Return a 1-D ndarray for vectorized single-component input."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'FDI']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 40.0]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))

    def test_array_values_match_scalar_calls(self):
        """Match vectorized values to equivalent scalar calls."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'FDI']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 40.0]),
        )
        expected_py = biomass.get_tree_biomass('PY', 1, 'wood', 30.0)
        expected_fdi = biomass.get_tree_biomass('FDI', 1, 'wood', 40.0)
        self.assertAlmostEqual(float(result[0]), expected_py, places=9)
        self.assertAlmostEqual(float(result[1]), expected_fdi, places=9)

    def test_array_multi_component_returns_tuple_of_ndarrays(self):
        """Return tuple[ndarray, ...] for vectorized multi-component requests."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'FDI']),
            np.array([1, 1]),
            ['wood', 'bark'],
            np.array([30.0, 40.0]),
        )
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], np.ndarray)
        self.assertIsInstance(result[1], np.ndarray)

    def test_mixed_scalar_array_input_returns_ndarray(self):
        """Broadcast scalar metadata across vector DBH input."""
        result = biomass.get_tree_biomass(
            'PY',
            1,
            'wood',
            np.array([30.0, 40.0]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))

    def test_array_with_height_uses_dbh_ht_coefficients(self):
        """Use DBH-height coefficients when ``height`` is provided."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'PY']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 30.0]),
            np.array([15.0, 15.0]),
        )
        py = biomass.BIOMASS_DATA['PY']
        expected = (
            biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['wood']]
            * py['BIOMASS_DBH-HT_Bwood1']
            * math.pow(30.0, py['BIOMASS_DBH-HT_Bwood2'])
            * math.pow(15.0, py['BIOMASS_DBH-HT_Bwood3'])
        )
        self.assertAlmostEqual(float(result[0]), expected, places=9)
        self.assertAlmostEqual(float(result[1]), expected, places=9)

    def test_array_mismatched_lengths_raises(self):
        """Raise ``ValueError`` for mismatched vector input lengths."""
        with self.assertRaises(ValueError):
            biomass.get_tree_biomass(
                np.array(['PY', 'FDI']),
                np.array([1]),
                'wood',
                np.array([30.0, 40.0]),
            )

    def test_array_no_height_uses_dbh_only_coefficients(self):
        """Use DBH-only coefficients for vector calls without ``height``."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'PY']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 40.0]),
        )
        py = biomass.BIOMASS_DATA['PY']
        decay_factor = biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['wood']]
        expected_0 = decay_factor * py['BIOMASS_DBH_Bwood1'] * math.pow(30.0, py['BIOMASS_DBH_Bwood2'])
        expected_1 = decay_factor * py['BIOMASS_DBH_Bwood1'] * math.pow(40.0, py['BIOMASS_DBH_Bwood2'])
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), expected_0, places=9)
        self.assertAlmostEqual(float(result[1]), expected_1, places=9)

    def test_hardwood_decay_class_out_of_range_raises(self):
        """Reject decay classes outside hardwood valid range."""
        with self.assertRaises(ValueError):
            biomass.get_tree_biomass('MA', 7, 'wood', 30.0)

    def test_vectorized_tree_biomass_returns_tuple_of_ndarrays(self):
        """Return expected vectorized outputs for mixed species/components."""
        result = biomass.get_tree_biomass(
            np.array(['PY', 'FDI']),
            np.array([1, 4]),
            ['wood', 'foliage'],
            np.array([30.0, 40.0]),
            np.array([15.0, 20.0]),
        )

        py = biomass.BIOMASS_DATA['PY']
        fdi = biomass.BIOMASS_DATA['FDI']
        expected_py_wood = (
            biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['wood']]
            * py['BIOMASS_DBH-HT_Bwood1']
            * math.pow(30.0, py['BIOMASS_DBH-HT_Bwood2'])
            * math.pow(15.0, py['BIOMASS_DBH-HT_Bwood3'])
        )
        expected_py_foliage = (
            biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['foliage']]
            * py['BIOMASS_DBH-HT_Bfoliage1']
            * math.pow(30.0, py['BIOMASS_DBH-HT_Bfoliage2'])
            * math.pow(15.0, py['BIOMASS_DBH-HT_Bfoliage3'])
        )
        expected_fdi_wood = (
            biomass.SOFTWOOD_DECAY[4][biomass.TREE_COMPONENT_INDEX['wood']]
            * fdi['BIOMASS_DBH-HT_Bwood1']
            * math.pow(40.0, fdi['BIOMASS_DBH-HT_Bwood2'])
            * math.pow(20.0, fdi['BIOMASS_DBH-HT_Bwood3'])
        )
        expected_fdi_foliage = (
            biomass.SOFTWOOD_DECAY[4][biomass.TREE_COMPONENT_INDEX['foliage']]
            * fdi['BIOMASS_DBH-HT_Bfoliage1']
            * math.pow(40.0, fdi['BIOMASS_DBH-HT_Bfoliage2'])
            * math.pow(20.0, fdi['BIOMASS_DBH-HT_Bfoliage3'])
        )

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], np.ndarray)
        self.assertIsInstance(result[1], np.ndarray)
        self.assertAlmostEqual(float(result[0][0]), expected_py_wood, places=9)
        self.assertAlmostEqual(float(result[1][0]), expected_py_foliage, places=9)
        self.assertAlmostEqual(float(result[0][1]), expected_fdi_wood, places=9)
        self.assertAlmostEqual(float(result[1][1]), expected_fdi_foliage, places=9)


class GetPhotoloadBiomassTests(unittest.TestCase):
    """Behavioral tests for ``biomass.get_photoload_biomass``."""

    def test_uses_default_height_when_missing(self):
        """Use species default height for scalar calls when ``height`` is omitted."""
        result = biomass.get_photoload_biomass('AMAL', 50.0)
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertAlmostEqual(result, expected, places=9)

    def test_scalar_input_returns_float(self):
        """Return a float for scalar photoload input."""
        result = biomass.get_photoload_biomass('AMAL', 50.0)
        self.assertIsInstance(result, float)
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertAlmostEqual(result, expected, places=9)

    def test_array_input_returns_ndarray(self):
        """Return a 1-D ndarray for vectorized photoload input."""
        result = biomass.get_photoload_biomass(
            np.array(['AMAL', 'VAGL']),
            np.array([50.0, 10.0]),
            np.array([35.56, 35.56]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))
        self.assertAlmostEqual(float(result[0]), 0.0148 * math.exp(0.0454 * 50.0), places=9)
        self.assertAlmostEqual(float(result[1]), 0.0052 * 10.0, places=9)

    def test_array_none_height_uses_all_defaults(self):
        """Apply default heights when ``height`` is omitted for arrays."""
        result = biomass.get_photoload_biomass(
            np.array(['AMAL', 'AMAL']),
            np.array([50.0, 50.0]),
        )
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), expected, places=9)
        self.assertAlmostEqual(float(result[1]), expected, places=9)

    def test_array_nan_height_uses_per_element_default(self):
        """Apply per-element default height where array height values are ``nan``."""
        result = biomass.get_photoload_biomass(
            np.array(['AMAL', 'AMAL']),
            np.array([50.0, 50.0]),
            np.array([np.nan, 20.0]),  # nan → default 35.56, explicit → 20.0
        )
        expected_default = 0.0148 * math.exp(0.0454 * 50.0)
        expected_explicit = (20.0 / 35.56) * 0.0148 * math.exp(0.0454 * 50.0)
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), expected_default, places=9)
        self.assertAlmostEqual(float(result[1]), expected_explicit, places=9)
        self.assertNotAlmostEqual(float(result[0]), float(result[1]), places=6)

    def test_invalid_code_in_array_warns_and_fills_zero(self):
        """Warn and emit zero biomass for invalid codes in vectorized calls."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = biomass.get_photoload_biomass(
                np.array(['AMAL', 'NOPE']),
                np.array([50.0, 20.0]),
            )
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), 0.0148 * math.exp(0.0454 * 50.0), places=9)
        self.assertAlmostEqual(float(result[1]), 0.0, places=9)
        self.assertTrue(caught)

    def test_array_mismatched_lengths_raises(self):
        """Raise ``ValueError`` for mismatched array lengths."""
        with self.assertRaises(ValueError):
            biomass.get_photoload_biomass(
                np.array(['AMAL', 'VAGL']),
                np.array([50.0]),
            )

    def test_invalid_photoload_species_warns_and_returns_zero(self):
        """Warn and return zero biomass for invalid scalar photoload species."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = biomass.get_photoload_biomass('NOPE', 20.0)

        self.assertEqual(result, 0.0)
        self.assertTrue(caught)


class GetDuffLitterBiomassTests(unittest.TestCase):
    """Behavioral tests for ``biomass.get_duff_litter_biomass``."""

    def test_bulk_density_for_single_species(self):
        """Return expected bulk density tuple for single-species input."""
        result = biomass.get_duff_litter_biomass('PY', return_type='bulk_density')
        self.assertEqual(result, (137.5986003, 53.50166766))

    def test_weighted_bulk_density_for_species_mix(self):
        """Return weighted bulk densities for mixed-species composition."""
        result = biomass.get_duff_litter_biomass(
            ['PY', 'FDI', 'PLI'],
            pct_list=[60.0, 30.0, 10.0],
            return_type='bulk_density',
        )
        expected_duff = (137.5986003 * 0.6) + (76.24788564 * 0.3) + (137.5986003 * 0.1)
        expected_litter = (53.50166766 * 0.6) + (29.63415723 * 0.3) + (53.50166766 * 0.1)
        self.assertAlmostEqual(result[0], expected_duff, places=9)
        self.assertAlmostEqual(result[1], expected_litter, places=9)

    def test_biomass_accepts_array_depths(self):
        """Support array depth inputs and return arrays for both layers."""
        result = biomass.get_duff_litter_biomass(
            'PY',
            return_type='biomass',
            duff_depth=np.array([3.4, 0.0]),
            litter_depth=np.array([1.4, 0.7]),
        )
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result[0], np.ndarray)
        self.assertIsInstance(result[1], np.ndarray)
        self.assertAlmostEqual(float(result[0][0]), 137.5986003 * 0.034, places=9)
        self.assertAlmostEqual(float(result[1][0]), 53.50166766 * 0.014, places=9)
        self.assertAlmostEqual(float(result[0][1]), 0.0, places=9)
        self.assertAlmostEqual(float(result[1][1]), 53.50166766 * 0.007, places=9)

    def test_scalar_biomass_returns_float(self):
        """Return float biomass for scalar depth input."""
        result = biomass.get_duff_litter_biomass('PY', return_type='biomass', duff_depth=3.4)
        self.assertIsInstance(result, float)
        self.assertAlmostEqual(result, 137.5986003 * 0.034, places=9)

    def test_array_both_depths_returns_tuple_of_ndarrays(self):
        """Return tuple of arrays when both duff and litter depths are arrays."""
        result = biomass.get_duff_litter_biomass(
            'PY',
            return_type='biomass',
            duff_depth=np.array([3.4, 5.0]),
            litter_depth=np.array([0.7, 1.0]),
        )
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result[0], np.ndarray)
        self.assertIsInstance(result[1], np.ndarray)
        self.assertAlmostEqual(float(result[0][0]), 137.5986003 * 0.034, places=9)
        self.assertAlmostEqual(float(result[1][0]), 53.50166766 * 0.007, places=9)
        self.assertAlmostEqual(float(result[0][1]), 137.5986003 * 0.050, places=9)
        self.assertAlmostEqual(float(result[1][1]), 53.50166766 * 0.010, places=9)

    def test_array_single_depth_returns_ndarray(self):
        """Return a single array when only one depth vector is provided."""
        result = biomass.get_duff_litter_biomass(
            'PY',
            return_type='biomass',
            duff_depth=np.array([3.4, 5.0]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), 137.5986003 * 0.034, places=9)
        self.assertAlmostEqual(float(result[1]), 137.5986003 * 0.050, places=9)

    def test_scalar_both_depths_returns_tuple_of_floats(self):
        """Return tuple of floats when both scalar depths are supplied."""
        result = biomass.get_duff_litter_biomass(
            'PY',
            return_type='biomass',
            duff_depth=3.4,
            litter_depth=0.7,
        )
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result[0], float)
        self.assertIsInstance(result[1], float)
        self.assertAlmostEqual(result[0], 137.5986003 * 0.034, places=9)
        self.assertAlmostEqual(result[1], 53.50166766 * 0.007, places=9)

    def test_requires_pct_list_for_species_mix(self):
        """Require ``pct_list`` when multiple species are provided."""
        with self.assertRaises(ValueError):
            biomass.get_duff_litter_biomass(['PY', 'FDI'], return_type='bulk_density')

    def test_mismatched_depth_array_lengths_raises(self):
        """Raise ``ValueError`` when depth vector lengths do not match."""
        with self.assertRaises(ValueError):
            biomass.get_duff_litter_biomass(
                'PY',
                return_type='biomass',
                duff_depth=np.array([3.4, 5.0]),
                litter_depth=np.array([0.7]),
            )

    def test_biomass_no_depth_raises(self):
        """Raise ``ValueError`` when biomass is requested without depths."""
        with self.assertRaises(ValueError):
            biomass.get_duff_litter_biomass('PY', return_type='biomass')


if __name__ == '__main__':
    unittest.main()
