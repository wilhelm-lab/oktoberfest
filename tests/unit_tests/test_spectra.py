import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from oktoberfest.data.spectra import Spectra, FragmentType


class TestSpectra(unittest.TestCase):
    """Class to test spectra."""

    def test_gen_column_names_pred(self):
        """Test gen_column_names pred."""
        columns = Spectra._gen_column_names(FragmentType.PRED)
        self.assertEqual(len(columns), 174)

    def test_gen_column_names_raw(self):
        """Test gen_column_names raw."""
        columns = Spectra._gen_column_names(FragmentType.RAW)
        self.assertEqual(len(columns), 174)

    def test_gen_column_names_mz(self):
        """Test gen_column_names mz."""
        columns = Spectra._gen_column_names(FragmentType.MZ)
        self.assertEqual(len(columns), 174)

    def test_resolve_prefix(self):
        """Test resolve_prefix."""
        self.assertEqual(Spectra._resolve_prefix(FragmentType.PRED), Spectra.INTENSITY_PRED_PREFIX)
        self.assertEqual(Spectra._resolve_prefix(FragmentType.RAW), Spectra.INTENSITY_COLUMN_PREFIX)
        self.assertEqual(Spectra._resolve_prefix(FragmentType.MZ), Spectra.MZ_COLUMN_PREFIX)

    def test_get_meta_data(self):
        """Test get_meta_data."""
        expected_result = pd.DataFrame({
            'metadata_1': ['a', 'b', 'c'],
            'metadata_2': ['d', 'e', 'f']
        })
        spectra = self.spectra_example()
        spectra.add_column(pd.Series(['d', 'e', 'f']), "metadata_2")
        result = spectra.get_meta_data()
        assert_frame_equal(expected_result, result)

    def test_add_matrix_from_hdf5(self):
        """Test add_matrix_from_hdf5."""
        intensity_data = np.random.randn(5, 174).tolist()
        spectra = Spectra()
        spectra.add_matrix_from_hdf5(pd.DataFrame(intensity_data), FragmentType.PRED)
        self.assertEqual(len(spectra.spectra_data.columns), 174)

    def spectra_example(self):
        """Example of spectra data."""
        spectra_data = pd.DataFrame({
            'INTENSITY_RAW_Y1+': [0.1, 0.2, 0.4],
            'INTENSITY_RAW_Y1++': [0.4, 0.5, 0.6],
            'MZ_RAW_Y1+': [100, 200, 300],
            'MZ_RAW_Y1++': [400, 500, 600],
            'INTENSITY_PRED_Y1+': [0.1, 0.2, 0.3],
            'INTENSITY_PRED_Y1++': [0.4, 0.5, 0.6],
            'metadata_1': ['a', 'b', 'c'],
        })
        spectra = Spectra()
        Spectra.add_columns(spectra, spectra_data)
        return spectra

