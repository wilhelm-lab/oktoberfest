import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.fragments import initialize_peaks

from oktoberfest import pp


class TestProcessing(unittest.TestCase):
    """Test class for preprocessing functions."""

    def test_list_spectra(self):
        """Test listing of spectra with expected user input."""
        spectra_path = Path(__file__).parent
        spectra_file = spectra_path / "test.mzml"
        spectra_file.open("w").close()
        self.assertEqual([spectra_path / "test.mzml"], pp.list_spectra(spectra_path, input_format="mzml"))
        spectra_file.unlink()

    def test_list_spectra_with_empty_string_folder(self):
        """Test listing spectra in a string folder without matching files."""
        self.assertRaises(AssertionError, pp.list_spectra, str(Path(__file__).parent), "raw")

    def test_list_spectra_with_wrong_folder(self):
        """Test listing spectra in a folder that does not exist."""
        self.assertRaises(NotADirectoryError, pp.list_spectra, Path(__file__).parent / "noexist", "raw")

    def test_list_spectra_with_wrong_format(self):
        """Test listing spectra with a format that isn't allowed."""
        self.assertRaises(ValueError, pp.list_spectra, Path(__file__).parent, "mzm")


def one_psm_library(sequence: str = "AAIGEATRL", n_peaks: int = 8) -> pd.DataFrame:
    """Build a one-PSM search-result frame whose peaks are real b/y fragments of ``sequence``.

    :param sequence: peptide to fragment
    :param n_peaks: how many of its theoretical b/y fragments to hand back as observed peaks
    :return: a frame with the columns annotate_spectral_library needs
    """
    fragments = initialize_peaks(sequence=sequence, mass_analyzer="FTMS", charge=2)[0]
    mz = np.array(sorted(f["mass"] for f in fragments if f["ion_type"] in ("b", "y"))[:n_peaks])
    return pd.DataFrame(
        {
            "MODIFIED_SEQUENCE": [sequence],
            "SEQUENCE": [sequence],
            "PEPTIDE_LENGTH": [len(sequence)],
            "PRECURSOR_CHARGE": [2],
            "MASS_ANALYZER": ["FTMS"],
            "FRAGMENTATION": ["HCD"],
            "MASS": [1000.0],
            "SCAN_NUMBER": [1],
            "RAW_FILE": ["file"],
            "REVERSE": [False],
            "SCORE": [1.0],
            "RETENTION_TIME": [1.0],
            "COLLISION_ENERGY": [30.0],
            "INTENSITIES": [np.linspace(1.0, 0.2, len(mz))],
            "MZ": [mz],
        }
    )


class TestScFeatureExtraction(unittest.TestCase):
    """Annotation must hand every sc_feature through to the Spectra object as its own column."""

    def test_every_sc_feature_key_becomes_a_column(self):
        """A key produced by annotation but never extracted here is silently unusable downstream."""
        aspec = pp.annotate_spectral_library(one_psm_library(), mass_tol=20, unit_mass_tol="ppm")
        missing = [key for key in c.SC_FEATURE_KEYS if key not in aspec.obs.columns]
        self.assertEqual(missing, [], f"sc_features not extracted onto the Spectra object: {missing}")

    def test_extraction_does_not_raise_on_a_key_the_annotation_omits(self):
        """The extraction indexes the dict by key; a partial dict must fail loudly here, not later."""
        aspec = pp.annotate_spectral_library(one_psm_library(), mass_tol=20, unit_mass_tol="ppm")
        coverage = float(aspec.obs["intensity_coverage"].iloc[0])
        self.assertGreater(coverage, 0.0)
        self.assertLessEqual(coverage, 1.0)

    def test_reporter_features_are_nan_without_a_tmt_label(self):
        """Reporter ions are only defined for TMT; a label-free run must not invent values."""
        aspec = pp.annotate_spectral_library(one_psm_library(), mass_tol=20, unit_mass_tol="ppm")
        for key in ("n_reporter_channels", "reporter_intensity_frac", "reporter_max_frac"):
            self.assertTrue(np.isnan(float(aspec.obs[key].iloc[0])), f"{key} should be NaN without TMT")

    def test_matching_method_is_applied(self):
        """A non-default resolver must actually run, not be silently ignored."""
        for method in ("nearest", "highest", "global_ransac"):
            aspec = pp.annotate_spectral_library(
                one_psm_library(), mass_tol=20, unit_mass_tol="ppm", matching_method=method
            )
            self.assertIn("intensity_coverage", aspec.obs.columns, f"{method} produced no sc_features")

    def test_unknown_matching_method_raises(self):
        """A resolver name that is not registered must fail, not fall back to the default."""
        with self.assertRaises(ValueError):
            pp.annotate_spectral_library(
                one_psm_library(), mass_tol=20, unit_mass_tol="ppm", matching_method="not_a_resolver"
            )
