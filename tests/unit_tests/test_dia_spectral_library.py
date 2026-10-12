import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c

from oktoberfest.data.spectra import FragmentType, Spectra
from oktoberfest.runner import _extract_gene_name, generate_spectral_lib_fdr_control
from oktoberfest.utils.config import Config


class TestExtractGeneName(unittest.TestCase):
    """Test gene name extraction from protein identifiers."""

    def test_extract_gene_name(self):
        """Test UniProt headers yield the gene name and other identifiers yield None."""
        self.assertEqual(_extract_gene_name("sp|P16402|H13_HUMAN"), "H13")
        self.assertEqual(_extract_gene_name("tr|A0A024R161|A0A024R161_HUMAN"), "A0A024R161")
        self.assertIsNone(_extract_gene_name("P16402"))
        self.assertIsNone(_extract_gene_name("sp|P16402"))
        self.assertIsNone(_extract_gene_name("sp|P16402|_HUMAN"))


class TestDIASpectralLibrary(unittest.TestCase):
    """Test the spectral library export after DIA rescoring."""

    def setUp(self):  # noqa: D102
        self.temp_dir = Path(tempfile.mkdtemp())
        self.data_dir = self.temp_dir / "out" / "data"
        self.fdr_dir = self.temp_dir / "out" / "results" / "percolator"
        self.data_dir.mkdir(parents=True)
        self.fdr_dir.mkdir(parents=True)

        self.config = Config()
        self.config.data["output"] = self.temp_dir / "out"
        self.config.data["fdr_estimation_method"] = "percolator"
        self.config.base_path = self.temp_dir
        # avoid validating the (not existing) search results and spectra inputs
        self.config.check = lambda: None

        self.var = Spectra._gen_vars_df()
        self.y1 = self.var.index.get_loc("y1+1")
        self.b2 = self.var.index.get_loc("b2+1")
        self.y3 = self.var.index.get_loc("y3+1")

    def tearDown(self):  # noqa: D102
        shutil.rmtree(self.temp_dir)

    def test_spectral_library_fdr_filter_and_columns(self):
        """Test only best PSMs per peptide and charge at 1% FDR are exported with their fragments."""
        self._write_pred_hdf5(
            "run_a",
            pd.DataFrame(
                {
                    "RAW_FILE": ["run_a"] * 4,
                    "SCAN_NUMBER": [1, 2, 3, 4],
                    "MODIFIED_SEQUENCE": ["PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "ACDEFK"],
                    "SEQUENCE": ["PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "ACDEFK"],
                    "PRECURSOR_CHARGE": [2, 2, 3, 2],
                    "MASS": [927.4549, 927.4549, 927.4549, 697.2737],
                    "PROTEINS": ["sp|P1|GENEA_HUMAN", "sp|P1|GENEA_HUMAN", "sp|P1|GENEA_HUMAN;sp|P2|GENEB_HUMAN", "P3"],
                    "RETENTION_TIME": [10.0, 12.0, 11.0, 20.0],
                    "retention_time_sec": [600.0, 720.0, 660.0, 1200.0],
                }
            ),
        )
        # sorted by score as written by percolator: the best PSM per peptide and charge is kept
        self._write_psms(
            [
                ("run_a-2-PEPTIDEK-2", "run_a", "_.PEPTIDEK._", 0.001),
                ("run_a-1-PEPTIDEK-2", "run_a", "_.PEPTIDEK._", 0.001),
                ("run_a-3-PEPTIDEK-3", "run_a", "_.PEPTIDEK._", 0.005),
                ("run_a-4-ACDEFK-2", "run_a", "_.ACDEFK._", 0.02),
            ]
        )

        generate_spectral_lib_fdr_control(self.config)

        lib = pd.read_csv(self.config.output / "results" / "spec_lib.tsv", sep="\t")
        precursors = lib[["ModifiedPeptideSequence", "PrecursorCharge", "NormalizedRetentionTime"]].drop_duplicates()
        # charge 2 from scan 2 (best PSM), charge 3 from scan 3, ACDEFK fails the FDR cutoff
        self.assertEqual(precursors.values.tolist(), [["PEPTIDEK", 2, 12.0], ["PEPTIDEK", 3, 11.0]])

        charge2 = lib[lib["PrecursorCharge"] == 2]
        # y3 (0.04) is below the intensity cutoff, b2 and y1 are kept, sorted by product m/z
        self.assertEqual(charge2["Annotation"].tolist(), ["y1+1", "b2+1"])
        np.testing.assert_allclose(charge2["LibraryIntensity"], [10000.0, 5000.0])
        np.testing.assert_allclose(charge2["ProductMz"], [147.11, 227.10])
        self.assertEqual(charge2["FragmentType"].tolist(), ["y", "b"])
        self.assertEqual(charge2["FragmentSeriesNumber"].tolist(), [1, 2])
        self.assertEqual(charge2["FragmentCharge"].tolist(), [1, 1])
        expected_mz = (927.4549 + 2 * c.PARTICLE_MASSES["PROTON"]) / 2
        np.testing.assert_allclose(charge2["PrecursorMz"], expected_mz)
        np.testing.assert_allclose(charge2["AverageExperimentalRetentionTime"], 660.0)
        self.assertTrue((charge2["ProteinId"] == "sp|P1|GENEA_HUMAN").all())
        self.assertTrue((charge2["GeneName"] == "GENEA").all())
        self.assertTrue((charge2["Proteotypic"] == 1).all())

        charge3 = lib[lib["PrecursorCharge"] == 3]
        self.assertTrue((charge3["AllMappedGenes"] == "GENEA;GENEB").all())
        self.assertTrue((charge3["Proteotypic"] == 0).all())
        np.testing.assert_allclose(charge3["PrecursorMz"], (927.4549 + 3 * c.PARTICLE_MASSES["PROTON"]) / 3)

    def test_spectral_library_with_scan_event_number(self):
        """Test charges are resolved correctly when the PSMId ends with the scan event number."""
        obs = pd.DataFrame(
            {
                "RAW_FILE": ["run_a"] * 2,
                "SCAN_NUMBER": [1, 2],
                "SCAN_EVENT_NUMBER": [5, 6],
                "MODIFIED_SEQUENCE": ["PEPTIDEK", "PEPTIDEK"],
                "SEQUENCE": ["PEPTIDEK", "PEPTIDEK"],
                "PRECURSOR_CHARGE": [2, 3],
                "MASS": [927.4549, 927.4549],
                "PROTEINS": ["P1", "P1"],
                "RETENTION_TIME": [10.0, 11.0],
                "retention_time_sec": [600.0, 660.0],
            }
        )
        self._write_pred_hdf5("run_a", obs)
        self._write_psms(
            [
                ("run_a-1-PEPTIDEK-2-5", "run_a", "_.PEPTIDEK._", 0.001),
                ("run_a-2-PEPTIDEK-3-6", "run_a", "_.PEPTIDEK._", 0.001),
            ]
        )

        generate_spectral_lib_fdr_control(self.config)

        lib = pd.read_csv(self.config.output / "results" / "spec_lib.tsv", sep="\t")
        self.assertEqual(sorted(lib["PrecursorCharge"].unique().tolist()), [2, 3])

    def test_spectral_library_requires_percolator(self):
        """Test the export fails for mokapot results."""
        self.config.data["fdr_estimation_method"] = "mokapot"
        with self.assertRaisesRegex(ValueError, "Percolator"):
            generate_spectral_lib_fdr_control(self.config)

    def _write_pred_hdf5(self, run_name: str, obs: pd.DataFrame):
        n = len(obs)
        mz = np.zeros((n, len(self.var)))
        mz[:, self.y1] = 147.11
        mz[:, self.b2] = 227.10
        mz[:, self.y3] = 375.20
        pred = np.zeros((n, len(self.var)))
        pred[:, self.y1] = 1.0
        pred[:, self.b2] = 0.5
        pred[:, self.y3] = 0.04
        spectra = Spectra(obs=obs, var=self.var)
        spectra.var_names = self.var.index
        spectra.add_mzs(mz, FragmentType.MZ)
        spectra.add_intensities_without_mapping(pred, FragmentType.PRED)
        spectra.write_as_hdf5(self.data_dir / f"{run_name}.mzml.pred.hdf5")

    def _write_psms(self, rows: list[tuple[str, str, str, float]]):
        pd.DataFrame(rows, columns=["PSMId", "filename", "peptide", "q-value"]).to_csv(
            self.fdr_dir / "rescore.percolator.psms.txt", sep="\t", index=False
        )
