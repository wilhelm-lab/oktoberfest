import unittest

import numpy as np
import pandas as pd

from oktoberfest import preprocessing as pp
from oktoberfest.data.spectra import FragmentType


def _dense(matrix) -> np.ndarray:
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


class TestAnnotateSpectralLibrary(unittest.TestCase):
    """Test the annotation of measured spectra into a Spectra object, without mocking the annotation step."""

    def test_annotate_linear(self):
        """Test the y1, b2 and y2 peaks of PEPTIDEK are annotated and stored under their ion names."""
        psms = pd.DataFrame(
            {
                "RAW_FILE": ["run1"],
                "SCAN_NUMBER": [1],
                "MODIFIED_SEQUENCE": ["PEPTIDEK"],
                "SEQUENCE": ["PEPTIDEK"],
                "PRECURSOR_CHARGE": [2],
                "PEPTIDE_LENGTH": [8],
                "MASS_ANALYZER": ["FTMS"],
                "FRAGMENTATION": ["HCD"],
                # y1+1, b2+1, y2+1 and one unmatched peak
                "INTENSITIES": [np.array([10.0, 20.0, 30.0, 5.0])],
                "MZ": [np.array([147.11280, 227.10263, 276.15539, 500.0])],
            }
        )
        library = pp.annotate_spectral_library(psms, mass_tol=20, unit_mass_tol="ppm")

        self.assertEqual(library.shape, (1, 174))
        intensities = pd.Series(_dense(library.get_matrix(FragmentType.RAW))[0], index=list(library.var_names))
        self.assertAlmostEqual(intensities["y2+1"], 1.0)
        self.assertAlmostEqual(intensities["b2+1"], 2 / 3)
        self.assertAlmostEqual(intensities["y1+1"], 1 / 3)
        self.assertIn("CALCULATED_MASS", library.obs.columns)

    def test_annotate_xl(self):
        """Test the y1 peaks of both peptides of a non-cleavable cross-link are annotated."""
        psms = pd.DataFrame(
            {
                "RAW_FILE": ["run1"],
                "SCAN_NUMBER": [1],
                "PRECURSOR_CHARGE": [4],
                "MASS_ANALYZER": ["FTMS"],
                "FRAGMENTATION": ["HCD"],
                "CROSSLINKER_TYPE": ["BS3"],
                "REVERSE": [False],
                "SEQUENCE_A": ["KSGKGFYIYQEGVK"],
                "SEQUENCE_B": ["VPPVQVSPLIK"],
                "CROSSLINKER_POSITION_A": [4],
                "CROSSLINKER_POSITION_B": [1],
                "PEPTIDE_LENGTH_A": [14],
                "PEPTIDE_LENGTH_B": [11],
                "MODIFIED_SEQUENCE_A": ["KSGK[UNIMOD:1898]GFYIYQEGVK"],
                "MODIFIED_SEQUENCE_B": ["K[UNIMOD:1898]PPVQVSPLIK"],
                # y1+1 of both peptides (K), y2+1 of peptide A (VK) and one unmatched peak
                "INTENSITIES": [np.array([10.0, 20.0, 30.0])],
                "MZ": [np.array([147.11280, 246.18121, 500.0])],
            }
        )
        library = pp.annotate_spectral_library_xl(psms, cms2=False, mass_tol=20, unit_mass_tol="ppm")

        self.assertEqual(library.shape, (1, 174))
        for fragment_type in (FragmentType.RAW_A, FragmentType.RAW_B):
            intensities = pd.Series(_dense(library.get_matrix(fragment_type))[0], index=list(library.var_names))
            self.assertGreater(intensities["y1+1"], 0.1, f"y1+1 not annotated for {fragment_type}")
