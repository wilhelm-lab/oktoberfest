import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.runner import _ce_calib, run_ce_calibration
from oktoberfest.utils import Config, ProcessStep


def _make_library(n: int = 2) -> Spectra:
    obs = pd.DataFrame({"SCAN_NUMBER": list(range(1, n + 1)), "COLLISION_ENERGY": [30.0] * n})
    obs.index = obs.index.astype(str)
    return Spectra(obs=obs, var=Spectra._gen_vars_df())


class SyncJobPool:
    """Stand-in for JobPool that runs jobs immediately in the current process."""

    def __init__(self, processes: int):  # noqa: D107
        self.processes = processes
        self.results = []

    def apply_async(self, func, args):  # noqa: D102
        self.results.append(func(*args))

    def check_pool(self):  # noqa: D102
        return self.results


@patch("oktoberfest.runner._get_best_ce")
@patch("oktoberfest.runner._annotate_and_get_library")
class TestCECalib(unittest.TestCase):
    """Test the per spectra file CE calibration step."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.config = Config()
        self.config.base_path = Path(self.tmp_dir.name)
        self.config.data["output"] = "out"
        self.config.data["inputs"] = {"spectra_type": "mzml"}
        (self.config.output / "data").mkdir(parents=True)
        (self.config.output / "proc").mkdir()
        self.spectra_file = Path(self.tmp_dir.name) / "run_a.mzML"
        self.pred_hdf5 = self.config.output / "data" / "run_a.mzml.pred.hdf5"

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def test_first_run_annotates_calibrates_and_stores_library(self, mock_annotate, mock_best_ce):
        """Test the library is annotated, calibrated, written and the step is marked as done."""
        library = _make_library()
        mock_annotate.return_value = library

        result = _ce_calib(self.spectra_file, self.config)

        self.assertIs(result, library)
        mock_annotate.assert_called_once_with(self.spectra_file, self.config, tims_meta_file=None)
        mock_best_ce.assert_called_once_with(library, self.spectra_file, self.config)
        self.assertTrue(self.pred_hdf5.is_file())
        self.assertTrue(ProcessStep(self.config.output, "ce_calib.run_a").is_done())

    def test_rerun_loads_stored_library(self, mock_annotate, mock_best_ce):
        """Test a finished calibration is not repeated and the stored library is returned."""
        _make_library(3).write_as_hdf5(self.pred_hdf5)
        ProcessStep(self.config.output, "ce_calib.run_a").mark_done()

        result = _ce_calib(self.spectra_file, self.config)

        self.assertEqual(result.n_obs, 3)
        mock_annotate.assert_not_called()
        mock_best_ce.assert_not_called()

    def test_rerun_with_missing_library_fails(self, mock_annotate, mock_best_ce):
        """Test a clear error if the step is marked as done but the stored library is missing."""
        ProcessStep(self.config.output, "ce_calib.run_a").mark_done()

        with self.assertRaisesRegex(FileNotFoundError, "run_a.mzml.pred.hdf5"):
            _ce_calib(self.spectra_file, self.config)
        mock_annotate.assert_not_called()

    def test_timstof_passes_metadata_file(self, mock_annotate, mock_best_ce):
        """Test timsTOF spectra are annotated with their timsmeta file."""
        mock_annotate.return_value = _make_library()
        for spectra_type, spectra_file in [("d", "run_a.d"), ("hdf", "run_a.hdf")]:
            with self.subTest(spectra_type=spectra_type):
                self.config.data["inputs"]["spectra_type"] = spectra_type
                (self.config.output / "proc" / "ce_calib.run_a.done").unlink(missing_ok=True)

                _ce_calib(Path(spectra_file), self.config)

                self.assertEqual(
                    mock_annotate.call_args.kwargs["tims_meta_file"], self.config.output / "msms" / "run_a.timsmeta"
                )


@patch("oktoberfest.runner._ce_calib")
@patch("oktoberfest.runner._preprocess")
@patch("oktoberfest.runner.pp.list_spectra")
@patch.object(Config, "check")
class TestRunCECalibration(unittest.TestCase):
    """Test the CollisionEnergyCalibration job."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.config = Config()
        self.config.base_path = Path(self.tmp_dir.name)
        self.config.data["output"] = "out"
        self.config.data["inputs"] = {"spectra": "spectra", "spectra_type": "mzml"}
        self.config.output.mkdir()
        self.spectra_files = [Path("run_a.mzML"), Path("run_b.mzML")]

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def test_calibrates_every_spectra_file(self, mock_check, mock_list_spectra, mock_preprocess, mock_ce_calib):
        """Test the config is validated, inputs are preprocessed and every spectra file is calibrated."""
        mock_list_spectra.return_value = self.spectra_files
        mock_preprocess.return_value = self.spectra_files

        run_ce_calibration(self.config)

        mock_check.assert_called_once()
        mock_preprocess.assert_called_once_with(self.spectra_files, self.config)
        self.assertEqual(
            mock_ce_calib.call_args_list,
            [call(self.spectra_files[0], self.config), call(self.spectra_files[1], self.config)],
        )

    @patch("oktoberfest.runner.JobPool", SyncJobPool)
    def test_calibrates_every_spectra_file_in_parallel(
        self, mock_check, mock_list_spectra, mock_preprocess, mock_ce_calib
    ):
        """Test every spectra file is calibrated when using multiple processes."""
        self.config.data["numThreads"] = 4
        mock_list_spectra.return_value = self.spectra_files
        mock_preprocess.return_value = self.spectra_files

        run_ce_calibration(self.config)

        self.assertEqual(sorted(c.args[0].name for c in mock_ce_calib.call_args_list), ["run_a.mzML", "run_b.mzML"])
