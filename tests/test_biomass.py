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
    def test_scalar_dbh_only_uses_dbh_coefficients(self):
        result = biomass.getTreeBiomass('PY', 1, 'wood', 30.0)
        expected = 0.0564 * math.pow(30.0, 2.4465)
        self.assertAlmostEqual(result, expected, places=9)

    def test_scalar_input_returns_float(self):
        result = biomass.getTreeBiomass('PY', 1, 'wood', 30.0)
        self.assertIsInstance(result, float)
        expected = 0.0564 * math.pow(30.0, 2.4465)
        self.assertAlmostEqual(result, expected, places=9)

    def test_array_single_component_returns_ndarray(self):
        result = biomass.getTreeBiomass(
            np.array(['PY', 'FDI']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 40.0]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))

    def test_array_values_match_scalar_calls(self):
        result = biomass.getTreeBiomass(
            np.array(['PY', 'FDI']),
            np.array([1, 1]),
            'wood',
            np.array([30.0, 40.0]),
        )
        expected_py = biomass.getTreeBiomass('PY', 1, 'wood', 30.0)
        expected_fdi = biomass.getTreeBiomass('FDI', 1, 'wood', 40.0)
        self.assertAlmostEqual(float(result[0]), expected_py, places=9)
        self.assertAlmostEqual(float(result[1]), expected_fdi, places=9)

    def test_array_multi_component_returns_tuple_of_ndarrays(self):
        result = biomass.getTreeBiomass(
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
        result = biomass.getTreeBiomass(
            'PY',
            1,
            'wood',
            np.array([30.0, 40.0]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))

    def test_array_with_height_uses_dbh_ht_coefficients(self):
        result = biomass.getTreeBiomass(
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
        with self.assertRaises(ValueError):
            biomass.getTreeBiomass(
                np.array(['PY', 'FDI']),
                np.array([1]),
                'wood',
                np.array([30.0, 40.0]),
            )

    def test_vectorized_tree_biomass_returns_tuple_of_ndarrays(self):
        result = biomass.getTreeBiomass(
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
    def test_uses_default_height_when_missing(self):
        result = biomass.getPhotoloadBiomass('AMAL', 50.0)
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertAlmostEqual(result, expected, places=9)

    def test_scalar_input_returns_float(self):
        result = biomass.getPhotoloadBiomass('AMAL', 50.0)
        self.assertIsInstance(result, float)
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertAlmostEqual(result, expected, places=9)

    def test_array_input_returns_ndarray(self):
        result = biomass.getPhotoloadBiomass(
            np.array(['AMAL', 'VAGL']),
            np.array([50.0, 10.0]),
            np.array([35.56, 35.56]),
        )
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, (2,))
        self.assertAlmostEqual(float(result[0]), 0.0148 * math.exp(0.0454 * 50.0), places=9)
        self.assertAlmostEqual(float(result[1]), 0.0052 * 10.0, places=9)

    def test_array_none_height_uses_all_defaults(self):
        result = biomass.getPhotoloadBiomass(
            np.array(['AMAL', 'AMAL']),
            np.array([50.0, 50.0]),
        )
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), expected, places=9)
        self.assertAlmostEqual(float(result[1]), expected, places=9)

    def test_array_nan_height_uses_per_element_default(self):
        result = biomass.getPhotoloadBiomass(
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
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = biomass.getPhotoloadBiomass(
                np.array(['AMAL', 'NOPE']),
                np.array([50.0, 20.0]),
            )
        self.assertIsInstance(result, np.ndarray)
        self.assertAlmostEqual(float(result[0]), 0.0148 * math.exp(0.0454 * 50.0), places=9)
        self.assertAlmostEqual(float(result[1]), 0.0, places=9)
        self.assertTrue(caught)

    def test_array_mismatched_lengths_raises(self):
        with self.assertRaises(ValueError):
            biomass.getPhotoloadBiomass(
                np.array(['AMAL', 'VAGL']),
                np.array([50.0]),
            )

    def test_invalid_photoload_species_warns_and_returns_zero(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = biomass.getPhotoloadBiomass('NOPE', 20.0)

        self.assertEqual(result, 0.0)
        self.assertTrue(caught)


class GetDuffLitterBiomassTests(unittest.TestCase):
    def test_bulk_density_for_single_species(self):
        result = biomass.getDuffLitterBiomass('PY', return_type='bulk_density')
        self.assertEqual(result, (137.5986003, 53.50166766))

    def test_weighted_bulk_density_for_species_mix(self):
        result = biomass.getDuffLitterBiomass(
            ['PY', 'FDI', 'PLI'],
            pct_list=[60.0, 30.0, 10.0],
            return_type='bulk_density',
        )
        expected_duff = (137.5986003 * 0.6) + (76.24788564 * 0.3) + (137.5986003 * 0.1)
        expected_litter = (53.50166766 * 0.6) + (29.63415723 * 0.3) + (53.50166766 * 0.1)
        self.assertAlmostEqual(result[0], expected_duff, places=9)
        self.assertAlmostEqual(result[1], expected_litter, places=9)

    def test_biomass_accepts_vector_depths(self):
        result = biomass.getDuffLitterBiomass(
            'PY',
            return_type='biomass',
            duff_depth=[3.4, 0.0],
            litter_depth=[0.7, 1.4],
        )
        expected = [
            (137.5986003 * 0.034, 53.50166766 * 0.007),
            (0.0, 53.50166766 * 0.014),
        ]
        self.assertEqual(len(result), 2)
        for actual_row, expected_row in zip(result, expected):
            self.assertAlmostEqual(actual_row[0], expected_row[0], places=9)
            self.assertAlmostEqual(actual_row[1], expected_row[1], places=9)

    def test_requires_pct_list_for_species_mix(self):
        with self.assertRaises(ValueError):
            biomass.getDuffLitterBiomass(['PY', 'FDI'], return_type='bulk_density')


if __name__ == '__main__':
    unittest.main()
