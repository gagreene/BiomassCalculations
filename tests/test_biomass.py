import math
import unittest
import warnings

import biomass


class GetTreeBiomassTests(unittest.TestCase):
    def test_scalar_dbh_only_uses_dbh_coefficients(self):
        result = biomass.getTreeBiomass('PY', 1, 'wood', 30.0)
        expected = 0.0564 * math.pow(30.0, 2.4465)
        self.assertAlmostEqual(result, expected, places=9)

    def test_vectorized_tree_biomass_returns_list_of_tuples(self):
        result = biomass.getTreeBiomass(
            ['PY', 'FDI'],
            [1, 4],
            ['wood', 'foliage'],
            [30.0, 40.0],
            [15.0, 20.0],
        )

        py = biomass.BIOMASS_DATA['PY']
        fdi = biomass.BIOMASS_DATA['FDI']
        expected = [
            (
                biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['wood']]
                * py['BIOMASS_DBH-HT_Bwood1']
                * math.pow(30.0, py['BIOMASS_DBH-HT_Bwood2'])
                * math.pow(15.0, py['BIOMASS_DBH-HT_Bwood3']),
                biomass.SOFTWOOD_DECAY[1][biomass.TREE_COMPONENT_INDEX['foliage']]
                * py['BIOMASS_DBH-HT_Bfoliage1']
                * math.pow(30.0, py['BIOMASS_DBH-HT_Bfoliage2'])
                * math.pow(15.0, py['BIOMASS_DBH-HT_Bfoliage3']),
            ),
            (
                biomass.SOFTWOOD_DECAY[4][biomass.TREE_COMPONENT_INDEX['wood']]
                * fdi['BIOMASS_DBH-HT_Bwood1']
                * math.pow(40.0, fdi['BIOMASS_DBH-HT_Bwood2'])
                * math.pow(20.0, fdi['BIOMASS_DBH-HT_Bwood3']),
                biomass.SOFTWOOD_DECAY[4][biomass.TREE_COMPONENT_INDEX['foliage']]
                * fdi['BIOMASS_DBH-HT_Bfoliage1']
                * math.pow(40.0, fdi['BIOMASS_DBH-HT_Bfoliage2'])
                * math.pow(20.0, fdi['BIOMASS_DBH-HT_Bfoliage3']),
            ),
        ]

        self.assertEqual(len(result), 2)
        for actual_row, expected_row in zip(result, expected):
            self.assertEqual(len(actual_row), 2)
            for actual, expected_value in zip(actual_row, expected_row):
                self.assertAlmostEqual(actual, expected_value, places=9)

    def test_rejects_mismatched_vector_lengths(self):
        with self.assertRaises(ValueError):
            biomass.getTreeBiomass(['PY', 'FDI'], [1], 'wood', [30.0, 40.0])


class GetPhotoloadBiomassTests(unittest.TestCase):
    def test_uses_default_height_when_missing(self):
        result = biomass.getPhotoloadBiomass('AMAL', 50.0)
        expected = 0.0148 * math.exp(0.0454 * 50.0)
        self.assertAlmostEqual(result, expected, places=9)

    def test_vectorized_photoload_biomass(self):
        result = biomass.getPhotoloadBiomass(['AMAL', 'VAGL'], [50.0, 10.0], [35.56, 35.56])
        expected = [
            0.0148 * math.exp(0.0454 * 50.0),
            0.0052 * 10.0,
        ]
        self.assertEqual(len(result), 2)
        for actual, expected_value in zip(result, expected):
            self.assertAlmostEqual(actual, expected_value, places=9)

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
