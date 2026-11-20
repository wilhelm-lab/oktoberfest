import inspect
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from oktoberfest.runner import input_xifdr, xl_psm_to_csm

from .. import DATA_DIR

XL_INPUTS = DATA_DIR / "xl" / "non-cleavable" / "inputs"
XIFDR_COLUMNS = [
    "run",
    "scan",
    "peptide1",
    "peptide2",
    "peptide link 1",
    "peptide link 2",
    "is decoy 1",
    "is decoy 2",
    "precursor charge",
    "accession1",
    "accession2",
    "peptide position 1",
    "peptide position 2",
    "score",
]


def _psm_to_csm(features_dir: str, original_or_rescore: str) -> pd.DataFrame:
    """Call xl_psm_to_csm for percolator results, with or without the percolator_or_mokapot argument."""
    if "percolator_or_mokapot" in inspect.signature(xl_psm_to_csm).parameters:
        return xl_psm_to_csm(features_dir, original_or_rescore, "percolator")
    return xl_psm_to_csm(features_dir, original_or_rescore)


def _spec_id(scan: int, decoy_1: bool, decoy_2: bool, index: int) -> str:
    """A cross-link SpecId as Oktoberfest writes it (the raw file name contains a dash, as in real data)."""
    return (
        f"run-A_1-{scan}-PEPK[UNIMOD:1898]TIDE-LK[UNIMOD:1898]R-3-decoy_p1-{decoy_1}-decoy_p2-{decoy_2}"
        f"-between-PEPKTIDE-LKR-{index}"
    )


class TestXlPsmToCsm(unittest.TestCase):
    """Test the conversion of cross-link percolator PSMs (one per peptide) to CSMs."""

    def test_rescore_percolator(self):
        """Test every CSM has both peptides, the lower score of the two, a TT/TD/DD label and no helper columns."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            for name in ("rescore.percolator.psms.txt", "rescore.percolator.decoy.psms.txt"):
                shutil.copy(XL_INPUTS / name, tmp_dir)
            csms = _psm_to_csm(tmp_dir, "rescore")

        psms = pd.concat(
            pd.read_csv(XL_INPUTS / name, sep="\t")
            for name in ("rescore.percolator.psms.txt", "rescore.percolator.decoy.psms.txt")
        )
        index = psms["PSMId"].str.rsplit("-", n=1).str[1].str.split("_")
        peptides_per_csm = index.str[0].groupby(index.str[0]).size()
        self.assertEqual(len(csms), int((peptides_per_csm == 2).sum()))
        self.assertTrue(csms["index_csm"].is_unique)
        self.assertTrue((csms["score"] == csms[["score_a", "score_b"]].min(axis=1)).all())
        self.assertTrue(set(csms["label"]) <= {"TT", "TD", "DD"})
        self.assertIn("TT", set(csms["label"]))
        self.assertIn("fdr_group", csms.columns)
        for helper in ("index", "_", "which_pep", "decoy_p1", "decoy_p2", "fdr_group_b"):
            self.assertNotIn(helper, csms.columns)


class TestInputXifdr(unittest.TestCase):
    """Test the xiFDR input written from percolator CSMs and the search engine's CSM features."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.fdr_dir = Path(self.tmp_dir.name)
        rows = [(101, False, False, 0.9, "p1_PROT1_p2_PROT2"), (102, False, False, 0.4, "p1_PROT3_p2_PROT1")]
        decoy_rows = [(103, True, False, -0.3, "p1_REV_PROT1_p2_PROT2")]
        for name, entries in (("rescore.percolator.csms.txt", rows), ("rescore.percolator.decoy.csms.txt", decoy_rows)):
            pd.DataFrame(
                {
                    "SpecId": [_spec_id(scan, d1, d2, i) for i, (scan, d1, d2, _, _) in enumerate(entries)],
                    "score": [score for *_, score, _ in entries],
                    "q-value": 0.01,
                    "Proteins": [proteins for *_, proteins in entries],
                }
            ).to_csv(self.fdr_dir / name, sep="\t", index=False)
        self.scans = [101, 102, 103]
        self.spec_ids = [f"run-A_1-{scan}-PEPK[UNIMOD:1898]TIDE-LK[UNIMOD:1898]R-3" for scan in self.scans]

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _write_features(self, columns: dict):
        pd.DataFrame({"SpecId": self.spec_ids, "SCORE": [10.0, 5.0, 1.0], **columns}).to_csv(
            self.fdr_dir / "rescore_features_csm.tab", sep="\t", index=False
        )

    def _check_percolator_input(self) -> pd.DataFrame:
        percolator = pd.read_csv(self.fdr_dir / "percolator_xifdr_input.csv").sort_values("scan").reset_index(drop=True)
        self.assertListEqual(list(percolator.columns), XIFDR_COLUMNS)
        self.assertListEqual(percolator["scan"].tolist(), self.scans)
        self.assertListEqual(percolator["score"].tolist(), [0.9, 0.4, -0.3])
        self.assertListEqual(percolator["accession1"].tolist(), ["PROT1", "PROT3", "REV_PROT1"])
        self.assertListEqual(percolator["accession2"].tolist(), ["PROT2", "PROT1", "PROT2"])
        self.assertListEqual(percolator["is decoy 1"].tolist(), [False, False, True])
        return percolator

    def test_xisearch(self):
        """Test the percolator and the xiSearch xiFDR inputs (percolator score vs. xiSearch match score)."""
        self._write_features(
            {
                "match_score": [10.0, 5.0, 1.0],
                "run_name": "run-A_1",
                "scan_number": self.scans,
                "sequence_p1": "PEPKTIDE",
                "sequence_p2": "LKR",
                "link_pos_p1": 4,
                "link_pos_p2": 2,
                "precursor_charge": 3,
                "start_pos_p1": [10, 20, 30],
                "start_pos_p2": [5, 6, 7],
            }
        )
        input_xifdr(str(self.fdr_dir), "xisearch")

        percolator = self._check_percolator_input()
        self.assertListEqual(percolator["peptide position 1"].tolist(), [10, 20, 30])
        xisearch = pd.read_csv(self.fdr_dir / "xisearch_xifdr_input.csv").sort_values("scan")
        self.assertListEqual(list(xisearch.columns), XIFDR_COLUMNS)
        self.assertListEqual(xisearch["score"].tolist(), [10.0, 5.0, 1.0])
        self.assertFalse((self.fdr_dir / "scout_xifdr_input.csv").exists())

    def test_scout(self):
        """Test the scout xiFDR input: peptides without the '+' of scout and the classification score."""
        self._write_features(
            {
                "ClassificationScore": [0.8, 0.6, 0.1],
                "FileName": "run-A_1",
                "ScanNumber": self.scans,
                "AlphaPeptide": "PEPK+TIDE",
                "BetaPeptide": "LK+R",
                "AlphaPos": 4,
                "BetaPos": 2,
                "Charge": 3,
                "Peptide Position 1": [10, 20, 30],
                "Peptide Position 2": [5, 6, 7],
            }
        )
        input_xifdr(str(self.fdr_dir), "scout")

        self._check_percolator_input()
        scout = pd.read_csv(self.fdr_dir / "scout_xifdr_input.csv").sort_values("scan")
        self.assertListEqual(list(scout.columns), XIFDR_COLUMNS)
        self.assertListEqual(scout["score"].tolist(), [0.8, 0.6, 0.1])
        self.assertListEqual(scout["peptide1"].tolist(), ["PEPKTIDE"] * 3)
        self.assertListEqual(scout["peptide2"].tolist(), ["LKR"] * 3)
