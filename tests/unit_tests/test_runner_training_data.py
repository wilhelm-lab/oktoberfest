import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from oktoberfest.data.spectra import FragmentType, Spectra
from oktoberfest.runner import _export_training_data_parquets
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

    def _write_spectra_hdf5(self, run_name: str, scan_numbers: list[int]):
        var_df = Spectra._gen_vars_df()
        obs = pd.DataFrame(
            {
                "RAW_FILE": [run_name] * len(scan_numbers),
                "SCAN_NUMBER": scan_numbers,
                "MODIFIED_SEQUENCE": [f"PEPTIDE{i}" for i in range(len(scan_numbers))],
                "PRECURSOR_CHARGE": [2] * len(scan_numbers),
                "COLLISION_ENERGY": [30.0] * len(scan_numbers),
                "FRAGMENTATION": ["HCD"] * len(scan_numbers),
            }
        )
        spectra = Spectra(obs=obs, var=var_df)
        spectra.var_names = var_df.index
        spectra.add_intensities_without_mapping(np.ones((len(scan_numbers), len(var_df))), FragmentType.RAW)
        spectra.write_as_hdf5(self.data_dir / f"{run_name}.mzml.hdf5")

    def _write_original_tab(self, rows: list[tuple[str, str]]):
        pd.DataFrame(rows, columns=["SpecId", "filename"]).to_csv(self.fdr_dir / "original.tab", sep="\t", index=False)

    def _write_original_psms(self, rows: list[tuple[str, float]]):
        pd.DataFrame(rows, columns=["PSMId", "q-value"]).to_csv(
            self.fdr_dir / "original.percolator.psms.txt", sep="\t", index=False
        )
