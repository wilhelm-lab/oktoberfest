import tempfile
import unittest
import unittest.mock
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c

from .. import DATA_DIR
from oktoberfest import rescore as re
from oktoberfest.data import Spectra

PREDICTION_FEATURES = [
    "spectral_angle",
    "abs_rt_diff",
    "count_predicted",
    "observed_intensity",
    "sum_observed_and_predicted",
]


class TestRescoring(unittest.TestCase):
    """Test class for rescoring functions."""

    def test_generate_features(self):
        """Test feature generation with internal results and additional columns."""
        cases = [
            ("library100.hdf5", "expected_rescore_all.tab", {"additional_columns": "all"}),
            ("library100.hdf5", "expected_rescore_none.tab", {"additional_columns": "none"}),
            ("library100.hdf5", "expected_rescore_list.tab", {"additional_columns": ["A"]}),
            (
                "library150_multifrag.hdf5",
                "expected_rescore_multifrag.tab",
                {"additional_columns": None, "task": "multifrag", "featured_ions": ["C", "z"]},
            ),
        ]
        for library_file, expected_file, kwargs in cases:
            with self.subTest(expected=expected_file), tempfile.TemporaryDirectory() as tmp_dir:
                # the retention time alignment model is cached next to the output file, keep it out of the test data
                output_file = Path(tmp_dir) / "rescore.tab"
                re.generate_features(
                    library=Spectra.from_hdf5(DATA_DIR / "library" / library_file),
                    search_type="rescore",
                    output_file=output_file,
                    all_features=False,
                    regression_method="spline",
                    **kwargs,
                )
                expected = pd.read_csv(DATA_DIR / "library" / expected_file, sep="\t")
                pd.testing.assert_frame_equal(expected, pd.read_csv(output_file, sep="\t"))

    def test_generate_features_original_and_rescore(self):
        """Test original features exclude prediction based features and both keep the percolator column layout."""
        library = Spectra.from_hdf5(DATA_DIR / "library" / "library100.hdf5")
        with tempfile.TemporaryDirectory() as tmp_dir:
            for search_type in ["original", "rescore"]:
                re.generate_features(
                    library=library, search_type=search_type, output_file=Path(tmp_dir) / f"{search_type}.tab"
                )
            original = pd.read_csv(Path(tmp_dir) / "original.tab", sep="\t")
            rescore = pd.read_csv(Path(tmp_dir) / "rescore.tab", sep="\t")

        for tab in (original, rescore):
            self.assertEqual(tab.columns[:5].tolist(), ["SpecId", "Label", "ScanNr", "filename", "ExpMass"])
            self.assertEqual(tab.columns[-2:].tolist(), ["Peptide", "Proteins"])
            self.assertEqual(len(tab), library.n_obs)
            self.assertFalse(tab.isna().any().any())
        self.assertTrue(set(PREDICTION_FEATURES).isdisjoint(original.columns))
        self.assertTrue(set(PREDICTION_FEATURES).issubset(rescore.columns))
        # the search engine score is a feature of both, its delta score only without predictions
        self.assertTrue({"andromeda", "andromeda_delta_score"}.issubset(original.columns))
        self.assertIn("andromeda", rescore.columns)
        self.assertNotIn("andromeda_delta_score", rescore.columns)
        shared = ["Charge2", "andromeda", "missedCleavages", "sequence_length"]
        pd.testing.assert_frame_equal(original[shared], rescore[shared])

        # ExpMass is the experimental precursor m/z, labels are 1 for targets and -1 for decoys
        obs = library.obs.reset_index(drop=True)
        expected_mz = (obs["MASS"] + obs["PRECURSOR_CHARGE"] * c.PARTICLE_MASSES["PROTON"]) / obs["PRECURSOR_CHARGE"]
        np.testing.assert_allclose(original["ExpMass"], expected_mz)
        np.testing.assert_array_equal(original["Label"], np.where(obs["REVERSE"], -1, 1))

    def test_generate_features_rejects_non_numeric_additional_column(self):
        """Test a clear error if a non-numeric search result column is requested as feature."""
        library = Spectra.from_hdf5(DATA_DIR / "library" / "library100.hdf5")
        with tempfile.TemporaryDirectory() as tmp_dir:
            with self.assertRaisesRegex(ValueError, "MASS_ANALYZER"):
                re.generate_features(
                    library=library,
                    search_type="original",
                    output_file=Path(tmp_dir) / "original.tab",
                    additional_columns=["MASS_ANALYZER"],
                )
            with self.assertRaisesRegex(ValueError, "NOT_A_COLUMN"):
                re.generate_features(
                    library=library,
                    search_type="original",
                    output_file=Path(tmp_dir) / "original.tab",
                    additional_columns=["NOT_A_COLUMN"],
                )


class TestMergeInput(unittest.TestCase):
    """Test merging the per spectra file tab files into the percolator / mokapot input."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp_dir.name)

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _write_tab_files(self) -> list[Path]:
        tabs = [
            pd.DataFrame({"SpecId": ["F1-1", "F1-2"], "Label": [1, -1], "spectral_angle": [0.7, np.nan]}),
            pd.DataFrame({"SpecId": ["F2-1"], "Label": [1], "spectral_angle": [0.5]}),
            pd.DataFrame({"SpecId": ["F3-1", "F3-2"], "Label": [-1, 1], "spectral_angle": [np.nan, 0.9]}),
        ]
        tab_files = []
        for i, tab in enumerate(tabs):
            tab_file = self.tmp_path / f"file{i}.rescore.tab"
            tab.to_csv(tab_file, sep="\t", index=False)
            tab_files.append(tab_file)
        return tab_files

    def test_merge_input(self):
        """Test files are concatenated in order with a single header and missing values set to 0."""
        for num_threads in [1, 2]:
            with self.subTest(num_threads=num_threads):
                tab_files = self._write_tab_files()
                output_file = self.tmp_path / f"merged_{num_threads}.tab"

                re.merge_input(tab_files=tab_files, output_file=output_file, num_threads=num_threads)

                merged = pd.read_csv(output_file, sep="\t")
                self.assertEqual(merged.columns.tolist(), ["SpecId", "Label", "spectral_angle"])
                self.assertEqual(merged["SpecId"].tolist(), ["F1-1", "F1-2", "F2-1", "F3-1", "F3-2"])
                self.assertEqual(merged["Label"].tolist(), [1, -1, 1, -1, 1])
                self.assertEqual(merged["spectral_angle"].tolist(), [0.7, 0.0, 0.5, 0.0, 0.9])
                # the per spectra file inputs are removed after merging
                self.assertFalse(any(tab_file.exists() for tab_file in tab_files))


class TestPercolatorFdrModified(unittest.TestCase):
    """Test the separate FDR estimation for unmodified and modified PSMs."""

    def test_split_and_ptm_features(self):
        """Test oxidation goes to the modified PSMs with a 0/1 feature per PTM, fixed modifications stay unmodified."""
        features = pd.DataFrame(
            {
                "SpecId": ["a", "b", "c", "d", "e"],
                "Label": [1, 1, -1, 1, -1],
                "spectral_angle": [0.9, 0.8, 0.3, 0.7, 0.2],
                "Peptide": [
                    "_.AC[UNIMOD:4]K._",
                    "_.AM[UNIMOD:35]K._",
                    "_.AS[UNIMOD:21]K._",
                    "_.AM[UNIMOD:35]S[UNIMOD:21]K._",
                    "_.R[UNIMOD:267]AK._",
                ],
                "Proteins": ["P"] * 5,
            }
        )
        with (
            tempfile.TemporaryDirectory() as tmp_dir,
            unittest.mock.patch.object(re.rescore, "rescore_with_percolator") as rescore,
        ):
            features_path = Path(tmp_dir) / "rescore.tab"
            features.to_csv(features_path, sep="\t", index=False)

            re.percolator_fdr_modified(features_path)

            unmodified = pd.read_csv(features_path, sep="\t")
            modified = pd.read_csv(Path(tmp_dir) / "mod_fdr" / "rescore.tab", sep="\t")
        self.assertEqual(rescore.call_count, 2)
        self.assertEqual(unmodified["SpecId"].tolist(), ["a", "e"])
        self.assertEqual(unmodified.columns.tolist(), features.columns.tolist())
        self.assertEqual(modified["SpecId"].tolist(), ["b", "c", "d"])
        self.assertEqual(
            modified.columns.tolist(),
            ["SpecId", "Label", "spectral_angle", "ptm_UNIMOD_21", "ptm_UNIMOD_35", "Peptide", "Proteins"],
        )
        self.assertEqual(modified["ptm_UNIMOD_21"].tolist(), [0, 1, 1])
        self.assertEqual(modified["ptm_UNIMOD_35"].tolist(), [1, 0, 1])

    def test_single_ptm_has_no_constant_feature(self):
        """Test a PTM carried by all modified PSMs adds no feature."""
        features = pd.DataFrame(
            {
                "SpecId": ["a", "b"],
                "Label": [1, -1],
                "spectral_angle": [0.9, 0.3],
                "Peptide": ["_.AM[UNIMOD:35]K._", "_.KM[UNIMOD:35]A._"],
                "Proteins": ["P", "P"],
            }
        )
        with (
            tempfile.TemporaryDirectory() as tmp_dir,
            unittest.mock.patch.object(re.rescore, "rescore_with_percolator"),
        ):
            features_path = Path(tmp_dir) / "rescore.tab"
            features.to_csv(features_path, sep="\t", index=False)

            re.percolator_fdr_modified(features_path)

            modified = pd.read_csv(Path(tmp_dir) / "mod_fdr" / "rescore.tab", sep="\t")
        self.assertEqual(modified.columns.tolist(), features.columns.tolist())
