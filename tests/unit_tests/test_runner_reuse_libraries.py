import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from oktoberfest.data import Spectra
from oktoberfest.runner import REUSED_CE_COLUMN, _library_from_previous_runs, _reuse_predictions
from oktoberfest.utils import Config

SPECTRA_FILE = Path("run1.mzML")


def _library(sequences: list[str], scans: list[int], ce: float, extra_column: str) -> Spectra:
    """Annotated and predicted library of a previous run: row i has raw intensity i + 1 and predicted intensity i + 0.5."""
    n = len(sequences)
    obs = pd.DataFrame(
        {
            "RAW_FILE": ["run1"] * n,
            "SCAN_NUMBER": scans,
            "MODIFIED_SEQUENCE": sequences,
            "PRECURSOR_CHARGE": [2] * n,
            "FRAGMENTATION": ["HCD"] * n,
            "COLLISION_ENERGY": [ce] * n,
            "PREDICTED_IRT": np.arange(n, dtype=float),
            extra_column: [1.0] * n,
        }
    )
    obs.index = obs.index.astype(str)
    library = Spectra(obs=obs, var=Spectra._gen_vars_df())
    values = np.zeros(library.shape)
    values[:, 0] = np.arange(n) + 1
    library.layers[Spectra.INTENSITY_LAYER_NAME] = csr_matrix(values)
    library.layers[Spectra.MZ_LAYER_NAME] = csr_matrix(values)
    values[:, 0] = np.arange(n) + 0.5
    library.layers[Spectra.INTENSITY_PRED_LAYER_NAME] = csr_matrix(values)
    return library


class TestReuseLibraries(unittest.TestCase):
    """Test reusing the annotation and predictions of previous runs (one per search engine) for merged results."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        root = Path(self.tmp_dir.name)
        self.runs = [root / "msfragger", root / "sage"]
        libraries = [
            _library(["AAAK", "CCCK"], [1, 2], ce=32.0, extra_column="EXPECT"),
            _library(["CCCK", "DDDK", "DDDK"], [2, 3, 3], ce=31.0, extra_column="SAGE_ONLY"),
        ]
        for run, library, best_ce in zip(self.runs, libraries, (32.0, 31.0), strict=True):
            (run / "data").mkdir(parents=True)
            (run / "results" / "ce_calibration").mkdir(parents=True)
            library.write_as_hdf5(run / "data" / "run1.mzml.pred.hdf5")
            with open(run / "results" / "ce_calibration" / "run1_ce.txt", "w") as f:
                json.dump({"27.0": best_ce}, f)
        self.config = Config()
        self.config.base_path = Path()
        self.config.data.update({"output": str(root / "merged"), "reuseLibraries": [str(run) for run in self.runs]})
        (root / "merged" / "msms").mkdir(parents=True)
        self.search = pd.DataFrame(
            {
                "RAW_FILE": ["run1"] * 3,
                "SCAN_NUMBER": [1, 2, 3],
                "MODIFIED_SEQUENCE": ["AAAK", "CCCK", "DDDK"],
                "PRECURSOR_CHARGE": [2] * 3,
            }
        )

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _library(self, search: pd.DataFrame | None = None) -> Spectra | None:
        with patch("oktoberfest.runner.pp.load_search", return_value=self.search if search is None else search):
            return _library_from_previous_runs(SPECTRA_FILE, self.config)

    def test_union_first_run_wins(self):
        """Test every candidate comes from the first run that has it, with the original CE and common columns."""
        library = self._library()
        self.assertEqual(list(library.obs["MODIFIED_SEQUENCE"]), ["AAAK", "CCCK", "DDDK"])
        self.assertEqual(list(library.obs[REUSED_CE_COLUMN]), [32.0, 32.0, 31.0])
        self.assertEqual(list(library.obs["COLLISION_ENERGY"]), [27.0, 27.0, 27.0])
        self.assertNotIn("EXPECT", library.obs)
        self.assertNotIn("SAGE_ONLY", library.obs)
        raw = library.layers[Spectra.INTENSITY_LAYER_NAME].toarray()[:, 0]
        np.testing.assert_array_equal(raw, [1.0, 2.0, 2.0])  # CCCK: msfragger row 2; DDDK: first of sage (row 2)

    def test_missing_candidate(self):
        """Test a candidate of this run missing in every previous run falls back to annotation (None)."""
        search = pd.concat([self.search, self.search.iloc[[0]].assign(SCAN_NUMBER=9)], ignore_index=True)
        self.assertIsNone(self._library(search))

    def test_missing_run_file(self):
        """Test a previous run without the library of the spectra file falls back to annotation (None)."""
        (self.runs[1] / "data" / "run1.mzml.pred.hdf5").unlink()
        self.assertIsNone(self._library())

    def test_reuse_predictions(self):
        """Test only candidates predicted at another CE are predicted again, in their rows."""
        library = self._library()
        library.obs["COLLISION_ENERGY"] = 32.0  # this run's calibration: the msfragger run's

        def predict(data, **kwargs):
            values = np.zeros(data.shape)
            values[:, 0] = 100.0
            data.layers[Spectra.INTENSITY_PRED_LAYER_NAME] = csr_matrix(values)

        with patch("oktoberfest.runner.pr.Predictor.from_config") as from_config:
            from_config.return_value.predict_intensities.side_effect = predict
            _reuse_predictions(library, SPECTRA_FILE, self.config)
        predicted = library.layers[Spectra.INTENSITY_PRED_LAYER_NAME].toarray()[:, 0]
        np.testing.assert_array_equal(predicted, [0.5, 1.5, 100.0])
        self.assertEqual(from_config.return_value.predict_intensities.call_args.kwargs["data"].n_obs, 1)
        self.assertNotIn(REUSED_CE_COLUMN, library.obs)
        np.testing.assert_array_equal(library.obs["PREDICTED_IRT"], [0.0, 1.0, 1.0])

    def test_ce_calibration_dir_default(self):
        """Test the CE calibration of the first previous run is used unless ce_calibration_dir is given."""
        self.assertEqual(self.config.ce_calibration_dir, self.runs[0] / "results" / "ce_calibration")
        self.config.data["ce_alignment_options"] = {"ce_calibration_dir": "other"}
        self.assertEqual(self.config.ce_calibration_dir, Path("other"))

    def test_config_check(self):
        """Test reuseLibraries must be a list of previous run folders."""
        self.config._check_reuse_libraries()
        self.config.data["reuseLibraries"] = "msfragger"
        with self.assertRaisesRegex(ValueError, "must be a list"):
            self.config._check_reuse_libraries()
        self.config.data["reuseLibraries"] = [str(self.runs[0] / "results")]
        with self.assertRaisesRegex(ValueError, "no data folder"):
            self.config._check_reuse_libraries()
