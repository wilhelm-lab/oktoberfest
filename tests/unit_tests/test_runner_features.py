import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.runner import _calculate_features
from oktoberfest.utils import Config, ProcessStep


def _make_library() -> Spectra:
    obs = pd.DataFrame({"SCAN_NUMBER": [1, 2], "PEPTIDE_LENGTH": [8, 9]})
    obs.index = obs.index.astype(str)
    return Spectra(obs=obs, var=Spectra._gen_vars_df())


@patch("oktoberfest.runner.re.generate_features")
@patch("oktoberfest.runner.pr.Predictor.from_config")
@patch("oktoberfest.runner._ce_calib")
class TestCalculateFeatures(unittest.TestCase):
    """Test prediction and feature generation per spectra file."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.config = Config()
        self.config.base_path = Path(self.tmp_dir.name)
        self.config.data["output"] = "out"
        self.config.data["models"] = {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"}
        self.config.data["fdr_estimation_method"] = "percolator"
        for sub_dir in ["data", "proc", "results"]:
            (self.config.output / sub_dir).mkdir(parents=True)
        self.spectra_file = Path("run_a.mzML")
        self.fdr_dir = self.config.output / "results" / "percolator"

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _written_tab_files(self, mock_generate_features) -> dict[str, Path]:
        return {c.kwargs["search_type"]: c.kwargs["output_file"] for c in mock_generate_features.call_args_list}

    def test_rescoring_predicts_and_writes_both_tab_files(self, mock_ce_calib, mock_from_config, mock_features):
        """Test predictions are made once and original and rescore features are generated."""
        library = _make_library()
        mock_ce_calib.return_value = library

        _calculate_features(self.spectra_file, self.config)

        self.assertEqual([c.kwargs["model_type"] for c in mock_from_config.call_args_list], ["intensity", "irt"])
        predictor = mock_from_config.return_value
        predictor.predict_intensities.assert_called_once()
        self.assertIsNone(predictor.predict_intensities.call_args.kwargs["chunk_idx"])
        predictor.predict_rt.assert_called_once_with(data=library)
        self.assertTrue((self.config.output / "data" / "run_a.mzml.pred.hdf5").is_file())
        self.assertEqual(
            self._written_tab_files(mock_features),
            {"original": self.fdr_dir / "run_a.original.tab", "rescore": self.fdr_dir / "run_a.rescore.tab"},
        )
        self.assertTrue(ProcessStep(self.config.output, "predict.run_a").is_done())
        self.assertTrue(ProcessStep(self.config.output, "calculate_features.run_a").is_done())

    def test_finished_predictions_are_not_repeated(self, mock_ce_calib, mock_from_config, mock_features):
        """Test predictions are skipped if the prediction step is done, but features are still generated."""
        mock_ce_calib.return_value = _make_library()
        ProcessStep(self.config.output, "predict.run_a").mark_done()

        _calculate_features(self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        self.assertEqual(set(self._written_tab_files(mock_features)), {"original", "rescore"})

    def test_finished_features_are_not_repeated(self, mock_ce_calib, mock_from_config, mock_features):
        """Test nothing is predicted or written if the feature step is done."""
        mock_ce_calib.return_value = _make_library()
        ProcessStep(self.config.output, "calculate_features.run_a").mark_done()

        _calculate_features(self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        mock_features.assert_not_called()

    def test_original_only_mode(self, mock_ce_calib, mock_from_config, mock_features):
        """Test GenerateTrainingData mode writes only original features without predictions in a separate step."""
        mock_ce_calib.return_value = _make_library()
        self.config._run_original = True

        _calculate_features(self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        self.assertEqual(self._written_tab_files(mock_features), {"original": self.fdr_dir / "run_a.original.tab"})
        self.assertTrue(ProcessStep(self.config.output, "calculate_features_original.run_a").is_done())
        # a later full rescoring in the same output folder is not skipped
        self.assertFalse(ProcessStep(self.config.output, "calculate_features.run_a").is_done())

    @patch("oktoberfest.runner.group_iterator")
    def test_alphapept_predictions_are_chunked_by_peptide_length(
        self, mock_group_iterator, mock_ce_calib, mock_from_config, mock_features
    ):
        """Test AlphaPept models get predictions chunked by peptide length."""
        library = _make_library()
        mock_ce_calib.return_value = library
        mock_group_iterator.return_value = iter([[0], [1]])
        self.config.data["models"]["intensity"] = "AlphaPept_ms2_generic"

        _calculate_features(self.spectra_file, self.config)

        self.assertEqual(mock_group_iterator.call_args.kwargs["group_by_column"], "PEPTIDE_LENGTH")
        self.assertEqual(mock_from_config.return_value.predict_intensities.call_args.kwargs["chunk_idx"], [[0], [1]])
