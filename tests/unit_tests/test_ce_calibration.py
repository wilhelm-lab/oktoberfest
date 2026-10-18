import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from oktoberfest.data.spectra import Spectra
from oktoberfest.runner import _get_best_ce
from oktoberfest.utils import Config

CE_RANGE = (25, 40)


def _make_library() -> Spectra:
    """
    Library with HCD spectra at several original CEs and CID spectra.

    The CID spectrum at CE 35 shares its CE with an HCD group, the CID spectrum at CE 33 shares its CE with the
    calibrated value of the HCD 30 group; neither may be touched by the HCD calibration.
    """
    obs = pd.DataFrame(
        {
            "FRAGMENTATION": ["HCD", "HCD", "HCD", "HCD", "CID", "CID"],
            "COLLISION_ENERGY": [30.0, 30.0, 33.0, 35.0, 35.0, 33.0],
            "PRECURSOR_CHARGE": [2, 3, 2, 3, 2, 3],
            "MASS": [1000.0, 1500.0, 1200.0, 1800.0, 1100.0, 1300.0],
        },
        index=[str(i) for i in range(6)],
    )
    return Spectra(obs=obs, var=Spectra._gen_vars_df())


def _make_alignment_library(group_column: str, best_ce_per_group: dict) -> SimpleNamespace:
    """Build an alignment library whose spectral angle peaks at the given best CE for each group."""
    rows = []
    for group_value, best_ce in best_ce_per_group.items():
        for ce in range(*CE_RANGE):
            rows.append(
                {group_column: group_value, "COLLISION_ENERGY": ce, "SPECTRAL_ANGLE": 1 - abs(ce - best_ce) / 100}
            )
    return SimpleNamespace(obs=pd.DataFrame(rows))


def _make_ransac_alignment_library(delta_ce: float) -> SimpleNamespace:
    """Build an alignment library in which the best CE is always the original CE shifted by delta_ce."""
    rows = []
    for orig_ce, charge, mass in [(30, 2, 1000.0), (30, 3, 1500.0), (33, 2, 1200.0), (35, 3, 1800.0), (33, 3, 1600.0)]:
        for ce in range(*CE_RANGE):
            rows.append(
                {
                    "ORIG_COLLISION_ENERGY": orig_ce,
                    "COLLISION_ENERGY": ce,
                    "PRECURSOR_CHARGE": charge,
                    "MASS": mass,
                    "SPECTRAL_ANGLE": 1 - abs(ce - (orig_ce + delta_ce)) / 100,
                }
            )
    return SimpleNamespace(obs=pd.DataFrame(rows))


@patch("oktoberfest.runner.pl")
@patch("oktoberfest.runner.pr.Predictor.from_config")
class TestCECalibration(unittest.TestCase):
    """Test class for testing ce calibration in the case of HCD, CID and ransac regression."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.config = Config()
        self.config.base_path = Path(self.tmp_dir.name)
        self.config.data["models"] = {"intensity": "Prosit_2020_intensity_HCD"}
        self.config.data["ce_alignment_options"] = {"ce_range": CE_RANGE}
        self.spectra_file = Path("file1.mzML")
        self.ce_file = self.config.output / "results" / "ce_calibration" / "file1_ce.txt"

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _set_mode(self, acquisition_method: str, use_ransac_model: bool):
        self.config.data["acquisitionMethod"] = acquisition_method
        self.config.data["ce_alignment_options"]["use_ransac_model"] = use_ransac_model

    def test_ce_calibration_cid(self, mock_from_config, mock_pl):
        """Test ce calibration bypassing for CID data."""
        self._set_mode("DDA", use_ransac_model=False)
        library = _make_library()
        library.obs["FRAGMENTATION"] = "CID"

        _get_best_ce(library, self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [35.0] * 6)
        self.assertEqual(json.loads(self.ce_file.read_text()), 35)

    def test_ce_calibration_only_calibrates_exact_hcd(self, mock_from_config, mock_pl):
        """Test spectra with other labels, e.g. EThcD or unknown, are not calibrated like HCD spectra."""
        self._set_mode("DDA", use_ransac_model=False)
        library = _make_library()
        library.obs["FRAGMENTATION"] = ["EThcD", "ETHCD", "unknown", "unknown", "CID", "CID"]

        _get_best_ce(library, self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [35.0] * 6)

    def test_ce_calibration_hcd(self, mock_from_config, mock_pl):
        """Test ce calibration by highest scoring target PSMS without ransac for HCD data (DDA, per original CE)."""
        self._set_mode("DDA", use_ransac_model=False)
        library = _make_library()
        # 30 -> 33 and 33 -> 35 make sure an assigned CE is not remapped again by another group
        best_ce_per_orig = {30.0: 33, 33.0: 35, 35.0: 31}
        predictor = mock_from_config.return_value
        predictor.ce_calibration.return_value = _make_alignment_library("ORIG_COLLISION_ENERGY", best_ce_per_orig)

        _get_best_ce(library, self.spectra_file, self.config)

        self.assertFalse(predictor.ce_calibration.call_args.args[2])  # group_by_charge
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [33.0, 33.0, 35.0, 31.0, 35.0, 35.0])
        self.assertEqual(json.loads(self.ce_file.read_text()), {"30.0": 33, "33.0": 35, "35.0": 31})

    def test_ce_calibration_reuses_previous_calibration(self, mock_from_config, mock_pl):
        """Test a previous calibration is applied without predictions if ce_calibration_dir is configured."""
        self._set_mode("DDA", use_ransac_model=False)
        previous_dir = Path(self.tmp_dir.name) / "previous"
        previous_dir.mkdir()
        (previous_dir / "file1_ce.txt").write_text(json.dumps({"30.0": 33, "33.0": 35, "35.0": 31}))
        self.config.data["ce_alignment_options"]["ce_calibration_dir"] = "previous"
        library = _make_library()

        _get_best_ce(library, self.spectra_file, self.config)

        mock_from_config.assert_not_called()
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [33.0, 33.0, 35.0, 31.0, 35.0, 35.0])
        self.assertEqual(json.loads(self.ce_file.read_text()), {"30.0": 33, "33.0": 35, "35.0": 31})

    def test_ce_calibration_missing_previous_calibration(self, mock_from_config, mock_pl):
        """Test the CE is calibrated as usual if the previous calibration has no file for this spectra file."""
        self._set_mode("DDA", use_ransac_model=False)
        (Path(self.tmp_dir.name) / "previous").mkdir()
        self.config.data["ce_alignment_options"]["ce_calibration_dir"] = "previous"
        library = _make_library()
        predictor = mock_from_config.return_value
        predictor.ce_calibration.return_value = _make_alignment_library("ORIG_COLLISION_ENERGY", {30.0: 31, 33.0: 34})

        _get_best_ce(library, self.spectra_file, self.config)

        predictor.ce_calibration.assert_called_once()
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [31.0, 31.0, 34.0, 35.0, 35.0, 35.0])

    def test_ce_calibration_hcd_uncalibrated_group(self, mock_from_config, mock_pl):
        """Test that HCD spectra whose original CE has no calibration keep their original CE."""
        self._set_mode("DDA", use_ransac_model=False)
        library = _make_library()
        predictor = mock_from_config.return_value
        predictor.ce_calibration.return_value = _make_alignment_library("ORIG_COLLISION_ENERGY", {30.0: 33})

        _get_best_ce(library, self.spectra_file, self.config)

        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [33.0, 33.0, 33.0, 35.0, 35.0, 35.0])

    def test_ce_calibration_dia(self, mock_from_config, mock_pl):
        """Test ce calibration for DIA data, which calibrates one CE per precursor charge."""
        self._set_mode("DIA", use_ransac_model=False)
        library = _make_library()
        predictor = mock_from_config.return_value
        predictor.ce_calibration.return_value = _make_alignment_library("PRECURSOR_CHARGE", {2: 27, 3: 29})

        _get_best_ce(library, self.spectra_file, self.config)

        self.assertTrue(predictor.ce_calibration.call_args.args[2])  # group_by_charge
        np.testing.assert_array_equal(library.obs["COLLISION_ENERGY"].to_numpy(), [27.0, 29.0, 27.0, 29.0, 35.0, 35.0])
        self.assertEqual(json.loads(self.ce_file.read_text()), {"2": 27, "3": 29})

    def test_ce_calibration_ransac(self, mock_from_config, mock_pl):
        """Test ce calibration with a ransac regressor for timsTOF data."""
        self._set_mode("DDA", use_ransac_model=True)
        library = _make_library()
        predictor = mock_from_config.return_value
        predictor.ce_calibration.return_value = _make_ransac_alignment_library(delta_ce=3)

        _get_best_ce(library, self.spectra_file, self.config)

        self.assertTrue(predictor.ce_calibration.call_args.args[2])  # group_by_charge
        # each HCD spectrum keeps its own original CE shifted by the regressed delta; CID spectra stay at 35
        np.testing.assert_allclose(
            library.obs["COLLISION_ENERGY"].to_numpy(), [33.0, 33.0, 36.0, 38.0, 35.0, 35.0], atol=1e-6
        )


if __name__ == "__main__":
    unittest.main()
