import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp
from spectrum_fundamentals.constants import EPSILON

from oktoberfest.rescore.chimeric import subtract_shared_peaks


class TestSubtractSharedPeaks(unittest.TestCase):
    """Test the removal of peaks of better candidates in chimeric spectra."""

    def setUp(self):
        """Two candidates of scan 1 share the peak at m/z 200 (the most intense one), scan 2 has one candidate."""
        self.meta = pd.DataFrame(
            {"RAW_FILE": ["f", "f", "f"], "SCAN_NUMBER": [1, 1, 2], "MOST_INTESE_PEAK": [1000.0, 2000.0, 500.0]}
        )
        # ions 0-3; candidate 0 matches 100 and 200 like its prediction, candidate 1 matches 200 and 300
        self.raw = np.array([[0.5, 1.0, EPSILON, EPSILON], [EPSILON, 1.0, 0.3, EPSILON], [1.0, 0.2, EPSILON, 0.0]])
        self.mz = np.array([[100.0, 200.0, 0, 0], [0, 200.0, 300.0, 0], [150.0, 250.0, 0, 0]])
        self.pred = np.array([[0.5, 1.0, EPSILON, EPSILON], [EPSILON, 0.2, 1.0, EPSILON], [1.0, 0.2, 0.1, 0.0]])

    def test_subtraction(self):
        """Test shared peaks are removed from worse candidates and their most intense remaining peak is used."""
        raw, mz, meta = subtract_shared_peaks(
            self.meta, sp.csr_matrix(self.raw), sp.csr_matrix(self.mz), sp.csr_matrix(self.pred)
        )
        raw, mz = raw.toarray(), mz.toarray()
        np.testing.assert_array_equal(raw[[0, 2]], self.raw[[0, 2]])
        np.testing.assert_array_equal(mz[[0, 2]], self.mz[[0, 2]])
        np.testing.assert_array_equal(raw[1], [EPSILON, EPSILON, 1.0, EPSILON])
        np.testing.assert_array_equal(mz[1], [0, 0, 300.0, 0])
        np.testing.assert_allclose(meta["MOST_INTESE_PEAK"], [1000.0, 600.0, 500.0])
        self.assertEqual(self.meta["MOST_INTESE_PEAK"].tolist(), [1000.0, 2000.0, 500.0])

    def test_most_intense_peak_kept(self):
        """Test intensities and the most intense peak stay as they are if the most intense peak is not shared."""
        raw = self.raw.copy()
        raw[1] = [EPSILON, 0.4, 1.0, EPSILON]
        out, _, meta = subtract_shared_peaks(
            self.meta, sp.csr_matrix(raw), sp.csr_matrix(self.mz), sp.csr_matrix(self.pred)
        )
        np.testing.assert_array_equal(out.toarray()[1], [EPSILON, EPSILON, 1.0, EPSILON])
        self.assertEqual(meta["MOST_INTESE_PEAK"].tolist(), [1000.0, 2000.0, 500.0])

    def test_single_candidates(self):
        """Test nothing changes without chimeric spectra."""
        meta = self.meta.assign(SCAN_NUMBER=[1, 2, 3])
        raw, mz, _ = subtract_shared_peaks(
            meta, sp.csr_matrix(self.raw), sp.csr_matrix(self.mz), sp.csr_matrix(self.pred)
        )
        np.testing.assert_array_equal(raw.toarray(), self.raw)
        np.testing.assert_array_equal(mz.toarray(), self.mz)
