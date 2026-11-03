import json
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from oktoberfest.runner import _speclib_from_digestion, prepare_rescore_xl_psm_level, xl_preprocessing_plot_csm
from oktoberfest.utils import Config


def _assert_no_chained_assignment(test: unittest.TestCase, caught: list[warnings.WarningMessage]):
    """Chained in-place assignments do nothing under Copy-on-Write (pandas 3) and warn under pandas 2."""
    messages = [str(w.message) for w in caught if "chained assignment" in str(w.message)]
    test.assertListEqual(messages, [])


class TestSpeclibFromPeptides(unittest.TestCase):
    """Test the spectral library generation from a peptide list (up to the internal format)."""

    def test_missing_proteins(self):
        """Test peptides without a protein get the protein "unknown"."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            (tmp_path / "peptides.csv").write_text("peptide,proteins\nPEPTIDEK,P1\nLGGNEQVTR,\nAGLLCPSTK,P2;P3\n")
            (tmp_path / "cfg_speclib.json").write_text(
                json.dumps(
                    {
                        "type": "SpectralLibraryGeneration",
                        "inputs": {"library_input": "peptides.csv", "library_input_type": "peptides"},
                        "output": "out",
                        "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
                        "spectralLibraryOptions": {"fragmentation": "HCD", "precursorCharge": [2], "nrOx": 0},
                    }
                )
            )
            config = Config()
            config.read(tmp_path / "cfg_speclib.json")
            config.output.mkdir()

            with (
                warnings.catch_warnings(record=True) as caught,
                patch("oktoberfest.runner.pp.gen_lib", side_effect=StopIteration) as mock_gen_lib,
            ):
                warnings.simplefilter("always")
                with self.assertRaises(StopIteration):
                    _speclib_from_digestion(config)
            _assert_no_chained_assignment(self, caught)

            mock_gen_lib.assert_called_once_with(config.output / "peptides_internal.csv")
            internal = pd.read_csv(config.output / "peptides_internal.csv")
            proteins = dict(zip(internal["modified_sequence"], internal["proteins"], strict=True))
            self.assertDictEqual(proteins, {"PEPTIDEK": "P1", "LGGNEQVTR": "unknown", "AGLLCPSTK": "P2;P3"})


class TestXLPreprocessingPlotCsm(unittest.TestCase):
    """Test the split of cross-link CSMs into targets (TT) and decoys (TD, DD)."""

    def test_target_decoy_split(self):
        """Test only TT CSMs are written as targets, TD and DD CSMs as decoys."""
        labels = ["TT", "TD", "DD", "TT"]
        df = pd.DataFrame(
            {
                "PSMId_a": [f"psm{i}" for i in range(len(labels))],
                "label": labels,
                "scan_number_a": range(len(labels)),
                "filename_a": "run1",
                "peptide_a": "PEPTIDEK",
                "score": [4.0, 3.0, 2.0, 1.0],
                "q-value": 0.01,
                "q-value_a": 0.01,
                "q-value_b": 0.01,
                "score_a": 1.0,
                "score_b": 1.0,
                "fdr_group_a": "self",
                "proteinIds_b": "P1",
            }
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            xl_preprocessing_plot_csm(tmp_dir, df, "rescore", "percolator")
            targets = pd.read_csv(Path(tmp_dir) / "rescore.percolator.csms.txt", sep="\t")
            decoys = pd.read_csv(Path(tmp_dir) / "rescore.percolator.decoy.csms.txt", sep="\t")

        self.assertListEqual(targets["SpecId"].tolist(), ["psm0", "psm3"])
        self.assertListEqual(decoys["SpecId"].tolist(), ["psm1", "psm2"])
        self.assertListEqual(targets["is_target"].tolist(), ["TT", "TT"])
        self.assertListEqual(decoys["is_target"].tolist(), ["TD", "DD"])


class TestPrepareRescoreXLPsmLevel(unittest.TestCase):
    """Test the conversion of a cross-link CSM tab file into one row per peptide."""

    def test_proteins(self):
        """
        Test the CSM rows are split into one row per peptide with proteins and target / decoy labels.

        Proteins are rebuilt from protein_p1 / protein_p2 after missing values are filled, so they are never missing.
        """
        tab = pd.DataFrame(
            {
                "SpecId": ["run1-1", "run1-2"],
                "Label": [1, -1],
                "ScanNr": [1, 2],
                "filename": "run1",
                "ExpMass": [1000.0, 1100.0],
                "Peptide": "PEPTIDEK",
                "Proteins": ["P1", None],
                "protein_p1": ["P1", None],
                "protein_p2": ["P2", "P3"],
                "decoy_p1": [False, True],
                "decoy_p2": [False, False],
                "fdr_group": ["between", "between"],
                "base_sequence_p1": ["PEPTIDEK", "KEDITPEP"],
                "base_sequence_p2": ["LGGNEQVTR", "LGGNEQVTR"],
                "spectral_angle_a": [0.9, 0.5],
                "spectral_angle_b": [0.8, 0.4],
            }
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            tab.to_csv(Path(tmp_dir) / "rescore.tab", sep="\t", index=False)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                psms = prepare_rescore_xl_psm_level(tmp_dir, "rescore")
        _assert_no_chained_assignment(self, caught)

        self.assertEqual(len(psms), 4)
        self.assertFalse(psms["Proteins"].isna().any())
        self.assertListEqual(psms["Label"].tolist(), [1, -1, 1, 1])
