import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from oktoberfest.data.spectra import FragmentType, Spectra
from oktoberfest.runner import _export_training_data_parquets, run_generate_training_data
from oktoberfest.utils.config import Config


class TestGenerateTrainingDataParquets(unittest.TestCase):
    """Test GenerateTrainingData parquet export helpers."""

    def setUp(self):  # noqa: D102
        self.temp_dir = Path(tempfile.mkdtemp())
        self.data_dir = self.temp_dir / "out" / "data"
        self.fdr_dir = self.temp_dir / "out" / "results" / "percolator"
        self.spectra_dir = self.temp_dir / "spectra"
        self.data_dir.mkdir(parents=True)
        self.fdr_dir.mkdir(parents=True)
        self.spectra_dir.mkdir()

        # mock config
        self.config = Config()
        self.config.data["output"] = self.temp_dir / "out"
        self.config.data["fdr_estimation_method"] = "percolator"
        self.config.base_path = self.temp_dir

    def tearDown(self):  # noqa: D102
        shutil.rmtree(self.temp_dir)

    def test_export_full_and_fdr1_parquets_for_multiple_runs(self):
        """Test full and 1% FDR parquet generation across multiple spectra files."""
        self._write_spectra_hdf5("run_a", [11, 12, 13])
        self._write_spectra_hdf5("run_b", [21, 22])
        self._write_original_tab(
            [
                ("run_a-11-PEPTIDEA-2-0", "run_a"),
                ("run_a-12-PEPTIDEB-2-1", "run_a"),
                ("run_a-13-PEPTIDEC-2-2", "run_a"),
                ("run_b-21-PEPTIDED-2-3", "run_b"),
                ("run_b-22-PEPTIDEE-2-4", "run_b"),
            ]
        )
        self._write_original_psms(
            [
                ("run_a-11-PEPTIDEA-2-0", 0.009),
                ("run_a-12-PEPTIDEB-2-1", 0.010),
                ("run_b-22-PEPTIDEE-2-4", 0.001),
            ]
        )

        _export_training_data_parquets(self.config)

        run_a_full = pd.read_parquet(self.data_dir / "run_a.parquet")
        run_a_fdr = pd.read_parquet(self.data_dir / "run_a.fdr1.parquet")
        run_b_full = pd.read_parquet(self.data_dir / "run_b.parquet")
        run_b_fdr = pd.read_parquet(self.data_dir / "run_b.fdr1.parquet")

        self.assertEqual(run_a_full["scan_number"].tolist(), [11, 12, 13])
        self.assertEqual(run_b_full["scan_number"].tolist(), [21, 22])
        self.assertEqual(run_a_fdr["scan_number"].tolist(), [11])
        self.assertEqual(run_b_fdr["scan_number"].tolist(), [22])

    def test_export_raises_for_original_tab_hdf5_row_mismatch(self):
        """Test clear failure when original.tab rows do not match the HDF5 rows for a run."""
        self._write_spectra_hdf5("run_a", [11, 12])
        self._write_original_tab([("run_a-11-PEPTIDEA-2-0", "run_a")])
        self._write_original_psms([("run_a-11-PEPTIDEA-2-0", 0.001)])

        with self.assertRaisesRegex(ValueError, "row count"):
            _export_training_data_parquets(self.config)

    def test_export_uses_original_tab_filenames_instead_of_spectra_directory(self):
        """Test exports are driven by original.tab run names rather than config spectra discovery."""
        (self.spectra_dir / "run_without_search_results.mzml").touch()
        self._write_spectra_hdf5("run_a", [11])
        self._write_original_tab([("run_a-11-PEPTIDEA-2-0", "run_a")])
        self._write_original_psms([("run_a-11-PEPTIDEA-2-0", 0.001)])

        _export_training_data_parquets(self.config)

        self.assertTrue((self.data_dir / "run_a.parquet").exists())
        self.assertFalse((self.data_dir / "run_without_search_results.parquet").exists())

    def test_export_parquet_columns_for_dlomix(self):
        """Test the parquet columns follow the format expected by DLOmix for training."""
        var_df = Spectra._gen_vars_df()
        # internal convention: -1 = fragment cannot exist, 0 = fragment not observed
        intensities = np.full((2, len(var_df)), -1.0)
        intensities[:, 0] = 1.0
        intensities[:, 1] = 0.0
        self._write_spectra_hdf5(
            "run_a",
            [11, 12],
            intensities=intensities,
            modified_sequences=["PEPTIDEK", "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"],
            charges=[2, 3],
            fragmentation=["HCD", "CID"],
        )
        self._write_original_tab([("run_a-11-PEPTIDEK-2-0", "run_a"), ("run_a-12-PEPTIDEK-3-1", "run_a")])
        self._write_original_psms([("run_a-11-PEPTIDEK-2-0", 0.001)])

        _export_training_data_parquets(self.config)

        full = pd.read_parquet(self.data_dir / "run_a.parquet")
        self.assertEqual(full["modified_sequence"].tolist(), ["[]-PEPTIDEK-[]", "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]-[]"])
        np.testing.assert_array_equal(
            np.stack(full["precursor_charge_onehot"]), [[0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]]
        )
        np.testing.assert_allclose(full["collision_energy_aligned_normed"], [0.3, 0.3])
        self.assertEqual(full["method_nbr"].tolist(), [2, 1])
        self.assertEqual(full["package"].tolist(), ["HCD", "CID"])
        self.assertEqual(full["raw_file"].tolist(), ["run_a", "run_a"])
        raw = np.stack(full["intensities_raw"])
        self.assertEqual(raw.shape, (2, len(var_df)))
        # observed intensities are kept, not observed fragments are 0 and impossible fragments -1
        np.testing.assert_array_equal(raw[:, 0], [1.0, 1.0])
        np.testing.assert_array_equal(raw[:, 1], [0.0, 0.0])
        self.assertTrue((raw[:, 2:] == -1).all())

    def test_run_generate_training_data_requires_percolator(self):
        """Test the job fails early for mokapot."""
        self.config.data["fdr_estimation_method"] = "mokapot"
        with patch.object(Config, "check"), patch("oktoberfest.runner.run_rescoring") as mock_rescoring:
            with self.assertRaisesRegex(ValueError, "percolator"):
                run_generate_training_data(self.config)
        mock_rescoring.assert_not_called()

    def test_run_generate_training_data_runs_original_rescoring_and_export(self):
        """Test rescoring is run with the original features only before exporting the parquet files."""
        run_original_during_rescoring = []
        with (
            patch.object(Config, "check"),
            patch(
                "oktoberfest.runner.run_rescoring",
                side_effect=lambda config: run_original_during_rescoring.append(config._run_original),
            ),
            patch("oktoberfest.runner._export_training_data_parquets") as mock_export,
        ):
            run_generate_training_data(self.config)

        self.assertEqual(run_original_during_rescoring, [True])
        mock_export.assert_called_once_with(config=self.config)
        self.assertFalse(self.config._run_original)

    def test_run_generate_training_data_resets_flag_on_failure(self):
        """Test the original-only flag is reset and nothing is exported if rescoring fails."""
        with (
            patch.object(Config, "check"),
            patch("oktoberfest.runner.run_rescoring", side_effect=RuntimeError("percolator failed")),
            patch("oktoberfest.runner._export_training_data_parquets") as mock_export,
        ):
            with self.assertRaises(RuntimeError):
                run_generate_training_data(self.config)

        mock_export.assert_not_called()
        self.assertFalse(self.config._run_original)

    def _write_spectra_hdf5(
        self,
        run_name: str,
        scan_numbers: list[int],
        intensities: np.ndarray = None,
        modified_sequences: list[str] = None,
        charges: list[int] = None,
        fragmentation: list[str] = None,
    ):
        var_df = Spectra._gen_vars_df()
        obs = pd.DataFrame(
            {
                "RAW_FILE": [run_name] * len(scan_numbers),
                "SCAN_NUMBER": scan_numbers,
                "MODIFIED_SEQUENCE": modified_sequences or [f"PEPTIDE{i}" for i in range(len(scan_numbers))],
                "PRECURSOR_CHARGE": charges or [2] * len(scan_numbers),
                "COLLISION_ENERGY": [30.0] * len(scan_numbers),
                "FRAGMENTATION": fragmentation or ["HCD"] * len(scan_numbers),
            }
        )
        if intensities is None:
            intensities = np.ones((len(scan_numbers), len(var_df)))
        spectra = Spectra(obs=obs, var=var_df)
        spectra.var_names = var_df.index
        spectra.add_intensities_without_mapping(intensities, FragmentType.RAW)
        spectra.write_as_hdf5(self.data_dir / f"{run_name}.mzml.hdf5")

    def _write_original_tab(self, rows: list[tuple[str, str]]):
        pd.DataFrame(rows, columns=["SpecId", "filename"]).to_csv(self.fdr_dir / "original.tab", sep="\t", index=False)

    def _write_original_psms(self, rows: list[tuple[str, float]]):
        pd.DataFrame(rows, columns=["PSMId", "q-value"]).to_csv(
            self.fdr_dir / "original.percolator.psms.txt", sep="\t", index=False
        )
