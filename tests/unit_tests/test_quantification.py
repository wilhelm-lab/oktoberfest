import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata
import numpy as np
import pandas as pd

from oktoberfest.utils import Config
from oktoberfest.utils import quantification as quant

RUNS = ["run_a", "run_b"]


def _ids() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "File Name": ["run_a", "run_b", "run_a", "run_b", "run_b"],
            "Full Sequence": ["PEPTIDEK/2", "PEPTIDEK/2", "AM[UNIMOD:35]K/3", "LLLK/2", "KKK/3"],
            "Protein Accession": ["P1_HUMAN", "P1_HUMAN", "P2_YEAST;P3_YEAST", "P4_ECOLI", "P5_ECOLI"],
        }
    )


def _peaks() -> pd.DataFrame:
    """One peak per reason: used, weaker peak, match-between-runs, shared peak, not detected and mass error."""
    return pd.DataFrame(
        {
            "File Name": ["run_a", "run_a", "run_b", "run_a", "run_b", "run_b", "run_a", "run_b"],
            "Full Sequence": [
                "PEPTIDEK/2",
                "PEPTIDEK/2",
                "PEPTIDEK/2",
                "AM[UNIMOD:35]K/3",
                "AM[UNIMOD:35]K/3",
                "LLLK/2|KKK/3",
                "LLLK/2",
                "AM[UNIMOD:35]K/3",
            ],  # fmt: skip
            "Peak intensity": [100.0, 40.0, 200.0, 50.0, 80.0, 30.0, 0.0, 70.0],
            "Peak Detection Type": ["MSMS", "MSMS", "MSMS", "MSMS", "MBR", "MSMS", "MSMS", "MSMS"],
            "Peak RT Start": [10.0, 30.0, 11.0, 20.0, 21.0, 40.0, np.nan, 22.0],
            "Peak RT Apex": [10.2, 30.2, 11.2, 20.2, 21.2, 40.2, np.nan, 22.2],
            "Peak RT End": [10.4, 30.4, 11.4, 20.4, 21.4, 40.4, np.nan, 22.4],
            "MS2 Retention Time": [10.1, 30.1, 11.1, 20.1, np.nan, 40.1, 50.0, 22.1],
            "Peak Apex Mass Error (ppm)": [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, np.nan, 12.0],
        }
    )


class TestHelpers(unittest.TestCase):
    """Test the small helpers of the MS1 quantification."""

    def test_target_proteins(self):
        """Test decoy proteins are removed from protein lists."""
        self.assertEqual(quant.target_proteins("P1;rev_P2;P3"), "P1;P3")
        self.assertEqual(quant.target_proteins("rev_P2"), "")

    def test_strip_modifications(self):
        """Test modifications, including terminal ones, are removed."""
        sequences = pd.Series(["AM[UNIMOD:35]K", "[UNIMOD:1]-AC[UNIMOD:4]K", "PEPTIDEK"])
        self.assertEqual(quant.strip_modifications(sequences).tolist(), ["AMK", "ACK", "PEPTIDEK"])

    def test_peptides_at_fdr(self):
        """Test the flanks of the Percolator peptides are removed and the q-value cutoff is applied."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "peptides.txt"
            path.write_text("PSMId\tq-value\tpeptide\na\t0.001\t_.AM[UNIMOD:35]K._\nb\t0.02\t_.PEPK._\n")
            self.assertEqual(quant.peptides_at_fdr(path, 0.01), {"AM[UNIMOD:35]K"})

    def test_split_protein_columns(self):
        """Test proteinIds are split into one protein per trailing column."""
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "psms.txt", Path(tmp) / "split.txt"
            src.write_text("PSMId\tscore\tproteinIds\nid1\t1.5\tP1;P2\nid2\t0.5\tP3\n")
            quant.split_protein_columns(src, dst)
            self.assertEqual(
                dst.read_text().splitlines(), ["PSMId\tscore\tproteinIds", "id1\t1.5\tP1\tP2", "id2\t0.5\tP3"]
            )


class TestFlashLFQInput(unittest.TestCase):
    """Test building the FlashLFQ identification table from the rescored PSMs."""

    def setUp(self):  # noqa: D102
        self.psms = pd.DataFrame(
            {
                "PSMId": ["run_a-10-AM[UNIMOD:35]K-2-5", "run_a-11-PEPK-3-6", "run_b-12-PEPK-3-7"],
                "filename": ["run_a", "run_a", "run_b"],
                "q-value": [0.001, 0.005, 0.2],
                "proteinIds": ["P1;rev_P9", "P2", "P2"],
            }
        )
        self.meta = pd.DataFrame(
            {
                "PSMId": ["run_a-10-AM[UNIMOD:35]K-2-5", "run_a-11-PEPK-3-6", "run_b-12-PEPK-3-7"],
                "MODIFIED_SEQUENCE": ["AM[UNIMOD:35]K", "PEPK", "PEPK"],
                "PRECURSOR_CHARGE": [2, 3, 3],
                "RETENTION_TIME": [10.5, 20.0, 21.0],
                "CALCULATED_MASS": [393.17, 471.27, 471.27],
            }
        )

    def test_input(self):
        """Test PSMs above the q-value cutoff are dropped and the full sequence is made charge specific."""
        ids = quant.flashlfq_input(self.psms, self.meta, 0.01)
        self.assertEqual(ids["Full Sequence"].tolist(), ["AM[UNIMOD:35]K/2", "PEPK/3"])
        self.assertEqual(ids["Base Sequence"].tolist(), ["AMK", "PEPK"])
        self.assertEqual(ids["Protein Accession"].tolist(), ["P1", "P2"])
        self.assertEqual(ids["Scan Retention Time"].tolist(), [10.5, 20.0])
        self.assertEqual(ids["Peptide Monoisotopic Mass"].tolist(), [393.17, 471.27])

    def test_missing_metadata(self):
        """Test PSMs without annotated spectra are reported."""
        with self.assertRaisesRegex(ValueError, "1 rescored PSMs"):
            quant.flashlfq_input(self.psms, self.meta.iloc[1:], 0.01)

    def test_psm_metadata(self):
        """Test the PSM ids are built as in the Percolator tables, with the scan event number if present."""
        obs = self.meta.drop(columns="PSMId").assign(
            RAW_FILE=["run_a", "run_a", "run_b"], SCAN_NUMBER=[10, 11, 12], SCAN_EVENT_NUMBER=[5, 6, 7]
        )
        obs.index = obs.index.astype(str)
        with tempfile.TemporaryDirectory() as tmp:
            anndata.AnnData(obs=obs).write_h5ad(Path(tmp) / "run_a.mzml.hdf5")
            meta = quant.psm_metadata(Path(tmp))
        self.assertEqual(meta["PSMId"].tolist(), self.meta["PSMId"].tolist())

    def test_psm_metadata_missing(self):
        """Test a clear error without annotated spectra."""
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(FileNotFoundError, "mzml.hdf5"):
                quant.psm_metadata(Path(tmp))


class TestPeakTable(unittest.TestCase):
    """Test which FlashLFQ peaks are used for the precursor intensities."""

    def test_reasons(self):
        """Test match-between-runs peaks are dropped and every other peak is used or has a reason."""
        rows = quant.peak_table(_peaks(), _ids(), ppm_filter=5)
        self.assertEqual(len(rows), 8)  # 7 identified peaks, the shared one listed for both precursors
        reasons = rows.set_index(["Modified sequence", "Charge", "Raw file", "Intensity"])["Reason"].to_dict()
        self.assertEqual(
            reasons,
            {
                ("AM[UNIMOD:35]K", 3, "run_a", 50.0): "",
                ("AM[UNIMOD:35]K", 3, "run_b", 70.0): "mass error",
                ("KKK", 3, "run_b", 30.0): "shared peak",
                ("LLLK", 2, "run_a", 0.0): "not detected",
                ("LLLK", 2, "run_b", 30.0): "shared peak",
                ("PEPTIDEK", 2, "run_a", 40.0): "weaker peak",
                ("PEPTIDEK", 2, "run_a", 100.0): "",
                ("PEPTIDEK", 2, "run_b", 200.0): "",
            },
        )
        self.assertEqual(rows["Used"].sum(), 3)
        self.assertEqual(rows.loc[rows["Intensity"] == 70.0, "Mass error to run median (ppm)"].item(), 10.0)
        self.assertEqual(rows.loc[rows["Intensity"] == 50.0, "Proteins"].item(), "P2_YEAST;P3_YEAST")

    def test_no_ppm_filter(self):
        """Test None keeps peaks with any mass error."""
        rows = quant.peak_table(_peaks(), _ids(), ppm_filter=None)
        self.assertTrue(rows.loc[rows["Intensity"] == 70.0, "Used"].item())


class TestNormalization(unittest.TestCase):
    """Test the normalization of the runs to the first one."""

    @staticmethod
    def _rows(ratios: list[float]) -> pd.DataFrame:
        n = len(ratios)
        return pd.DataFrame(
            {
                "Raw file": ["run_a"] * n + ["run_b"] * n,
                "Modified sequence": [f"PEP{i}K" for i in range(n)] * 2,
                "Charge": 2,
                "Intensity": [1000.0] * n + [1000.0 * 2**r for r in ratios],
                "Used": True,
            }
        )

    def test_mode_and_median(self):
        """Test the mode follows the unchanged majority and the median is pulled by the changed precursors."""
        rows = self._rows([0.5] * 60 + [-2.0] * 40)
        mode = quant.normalization_offsets(rows, RUNS, "mode")
        median = quant.normalization_offsets(rows, RUNS, "median")
        self.assertEqual(mode["run_a"], 0.0)
        self.assertAlmostEqual(mode["run_b"], 0.5)
        self.assertAlmostEqual(median["run_b"], 0.5)
        rows = self._rows([0.5] * 40 + [-2.0] * 30 + [1.0] * 30)
        self.assertAlmostEqual(quant.normalization_offsets(rows, RUNS, "mode")["run_b"], 0.5)
        self.assertAlmostEqual(quant.normalization_offsets(rows, RUNS, "median")["run_b"], 0.5)
        # more than half of the precursors change (spread out): the median moves, the mode stays at 0.5
        rows = self._rows([0.5] * 40 + np.linspace(-2.5, -0.5, 55).tolist() + [1.5] * 5)
        self.assertAlmostEqual(quant.normalization_offsets(rows, RUNS, "mode")["run_b"], 0.5)
        self.assertLess(quant.normalization_offsets(rows, RUNS, "median")["run_b"], -0.5)
        self.assertEqual(quant.normalization_offsets(rows, RUNS, "none").tolist(), [0.0, 0.0])

    def test_no_shared_precursors(self):
        """Test a run without precursors in common with the reference is not shifted."""
        rows = quant.peak_table(_peaks(), _ids(), ppm_filter=5)
        rows = rows[~((rows["Raw file"] == "run_b") & (rows["Modified sequence"] == "PEPTIDEK"))]
        self.assertEqual(quant.normalization_offsets(rows, RUNS, "mode").tolist(), [0.0, 0.0])

    def test_normalization_table(self):
        """Test the summary lists every run with its offset and the precursors shared with the reference."""
        rows = quant.peak_table(_peaks(), _ids(), ppm_filter=5)
        table = quant.normalization_table(rows, RUNS, pd.Series({"run_a": 0.0, "run_b": 1.0}), "median")
        self.assertEqual(table["Raw file"].tolist(), RUNS)
        self.assertEqual(table["Log2 offset"].tolist(), [0.0, 1.0])
        self.assertEqual(table["Precursors shared with reference"].tolist(), [2, 1])


class TestDriftNormalization(unittest.TestCase):
    """Test the correction of the intensity drift along the retention time."""

    @staticmethod
    def _rows() -> pd.DataFrame:
        """400 precursors in 3 runs; run_c is 1 log2 higher after 50 min and has 30% precursors changed by +2."""
        n = 400
        rt = np.linspace(1.0, 99.0, n)
        base = 1000.0 * 2 ** np.sin(np.arange(n))
        log2_change = np.where(rt > 50, 1.0, 0.0) + np.where(np.arange(n) % 10 < 3, 2.0, 0.0)
        return pd.DataFrame(
            {
                "Raw file": np.repeat(["run_a", "run_b", "run_c"], n),
                "Modified sequence": [f"PEP{i}K" for i in range(n)] * 3,
                "Charge": 2,
                "Intensity": np.concatenate([base, base, base * 2**log2_change]),
                "Peak RT apex": np.tile(rt, 3),
                "Used": True,
            }
        )

    def test_curves(self):
        """Test the curve follows the drift of the unchanged majority and is flat for the runs without drift."""
        curves = quant.drift_curves(self._rows(), ["run_a", "run_b", "run_c"], min_per_bin=50)
        self.assertEqual(
            curves.columns.tolist(), ["Raw file", "Normalization", "Retention time", "Log2 offset", "Precursors"]
        )
        self.assertTrue((curves["Precursors"] >= 50).all())
        self.assertTrue(np.allclose(curves.loc[curves["Raw file"] != "run_c", "Log2 offset"], 0.0))
        run_c = curves[curves["Raw file"] == "run_c"]
        self.assertTrue(np.allclose(run_c.loc[run_c["Retention time"] < 30, "Log2 offset"], 0.0))
        self.assertTrue(np.allclose(run_c.loc[run_c["Retention time"] > 70, "Log2 offset"], 1.0))

    def test_correction(self):
        """Test the corrected intensities of the unchanged precursors agree between the runs."""
        rows = self._rows()
        curves = quant.drift_curves(rows, ["run_a", "run_b", "run_c"], min_per_bin=50)
        corrected = rows["Intensity"] / np.power(2.0, quant.drift_correction(rows, curves))
        unchanged = (np.arange(len(rows)) % 400 % 10 >= 3) & ((rows["Peak RT apex"] < 30) | (rows["Peak RT apex"] > 70))
        log = np.log2(corrected[unchanged]).to_numpy().reshape(3, -1)
        self.assertTrue(np.allclose(log[2], log[0]))

    def test_correction_without_curve_or_apex(self):
        """Test peaks without an apex or without a curve for their run are not corrected."""
        rows = pd.DataFrame({"Raw file": ["run_a", "run_a", "run_b"], "Peak RT apex": [5.0, np.nan, 5.0]})
        curves = pd.DataFrame({"Raw file": "run_a", "Retention time": [0.0, 10.0], "Log2 offset": [0.0, 1.0]})
        self.assertEqual(quant.drift_correction(rows, curves).tolist(), [0.5, 0.0, 0.0])

    def test_merge_bins(self):
        """Test bins are merged until they are large enough and a short last bin joins the previous one."""
        merged = quant._merge_bins(np.array([0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 3.0]), 3)
        self.assertEqual(merged.tolist(), [0, 0, 0, 1, 1, 1, 1])


class TestPrecursorTable(unittest.TestCase):
    """Test the precursor and peptidoform tables."""

    def setUp(self):  # noqa: D102
        self.rows = quant.peak_table(_peaks(), _ids(), ppm_filter=5)
        offsets = quant.normalization_offsets(self.rows, RUNS, "median")
        self.rows["Normalized intensity"] = self.rows["Intensity"] / np.power(2.0, self.rows["Raw file"].map(offsets))

    def test_precursors(self):
        """Test one row per precursor with the normalized intensity of its used peak in every run."""
        table = quant.precursor_table(self.rows, RUNS)
        self.assertEqual(table.columns.tolist(), ["Sequence", "Modified sequence", "Proteins", "Charge"] + RUNS)
        self.assertEqual(table["Modified sequence"].tolist(), ["AM[UNIMOD:35]K", "PEPTIDEK"])
        self.assertEqual(table["Sequence"].tolist(), ["AMK", "PEPTIDEK"])
        self.assertEqual(table["Charge"].tolist(), [3, 2])
        self.assertEqual(table["Proteins"].tolist(), ["P2_YEAST;P3_YEAST", "P1_HUMAN"])
        self.assertEqual(table["run_a"].tolist(), [50.0, 100.0])
        self.assertTrue(np.isnan(table["run_b"].iloc[0]))  # MBR and mass error peaks are not used
        self.assertEqual(table["run_b"].iloc[1], 100.0)  # 200 shifted by the median log2 ratio of 1

    def test_peptidoforms(self):
        """Test the charge states of a modified sequence are summed, and a run without any stays missing."""
        precursors = pd.DataFrame(
            {
                "Sequence": ["PEPK", "PEPK", "AMK"],
                "Modified sequence": ["PEPK", "PEPK", "AM[UNIMOD:35]K"],
                "Proteins": ["P1", "P1", "P2"],
                "Charge": [2, 3, 2],
                "run_a": [10.0, 5.0, np.nan],
                "run_b": [np.nan, 7.0, 3.0],
            }
        )
        table = quant.peptidoform_table(precursors)
        self.assertEqual(table.columns.tolist(), ["Sequence", "Modified sequence", "Proteins"] + RUNS)
        self.assertEqual(table["Modified sequence"].tolist(), ["AM[UNIMOD:35]K", "PEPK"])
        self.assertTrue(np.isnan(table["run_a"].iloc[0]))
        self.assertEqual(table["run_a"].iloc[1], 15.0)
        self.assertEqual(table["run_b"].tolist(), [3.0, 7.0])
        filtered = quant.peptidoform_table(precursors, peptides={"PEPK"})
        self.assertEqual(filtered["Modified sequence"].tolist(), ["PEPK"])

    def test_precursors_by_protein_group(self):
        """Test precursors are assigned to their unique accepted protein group, others are dropped."""
        rows = self.rows.assign(**{"Normalized intensity": self.rows["Intensity"]})
        table = quant.precursor_table(rows, RUNS)
        groups = pd.DataFrame({"Protein IDs": ["P1_HUMAN", "P2_YEAST", "P3_YEAST"]})
        aq = quant.precursors_by_protein_group(table, groups)
        self.assertEqual(aq.columns.tolist(), ["protein", "ion"] + RUNS)
        self.assertEqual(aq["ion"].tolist(), ["PEPTIDEK/2"])  # AMK maps to two groups
        self.assertEqual(aq[RUNS].values.tolist(), [[100.0, 200.0]])

        groups = pd.DataFrame({"Protein IDs": ["P1_HUMAN", "P2_YEAST;P3_YEAST"]})
        aq = quant.precursors_by_protein_group(table, groups)
        self.assertEqual(aq["protein"].tolist(), ["P2_YEAST;P3_YEAST", "P1_HUMAN"])
        self.assertEqual(aq["run_b"].tolist(), [0.0, 200.0])


class TestRunFlashLFQ(unittest.TestCase):
    """Test the FlashLFQ call."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp_dir.name)
        self.config = Config()
        self.config.base_path = self.dir
        self.config.data["numThreads"] = 3
        self.config.data["quantificationOptions"] = {"flashlfqExe": "CMD.dll", "dotnetExe": "/opt/dotnet"}

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    @patch("oktoberfest.utils.quantification.subprocess.run")
    def test_command(self, mock_run):
        """Test the .dll is run with dotnet, with the identified charge only, and a missing output is an error."""
        mock_run.return_value = MagicMock(returncode=0)
        with self.assertRaisesRegex(RuntimeError, "FlashLFQ failed"):
            quant.run_flashlfq(self.config, self.dir)
        cmd = mock_run.call_args.args[0]
        self.assertEqual(cmd[:2], ["/opt/dotnet", str(self.dir / "CMD.dll")])
        self.assertIn("--chg", cmd)
        self.assertEqual(cmd[cmd.index("--thr") + 1], "3")
        self.assertNotIn("--nor", cmd)

        mock_run.side_effect = lambda *args, **kwargs: (
            (self.dir / "QuantifiedPeaks.tsv").touch() or MagicMock(returncode=0)
        )
        quant.run_flashlfq(self.config, self.dir)


class TestApplyQuant(unittest.TestCase):
    """Test the quantification route is chosen by quantificationOptions.precursor."""

    @patch("oktoberfest.utils.quantification.quantify_search_engine")
    @patch("oktoberfest.utils.quantification.quantify_ms1")
    def test_dispatch(self, mock_ms1, mock_search_engine):
        """Test FlashLFQ runs the MS1 quantification and the default the search engine's quantification."""
        config = Config()
        quant.apply_quant(config)
        mock_search_engine.assert_called_once_with(config)
        config.data["quantificationOptions"] = {"precursor": "flashlfq"}
        quant.apply_quant(config)
        mock_ms1.assert_called_once_with(config)
