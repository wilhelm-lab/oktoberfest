import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from oktoberfest.rescore import add_adjacent_decoy_candidates, localization_groups, neighbor_classes


class TestAdjacentDecoys(unittest.TestCase):
    """Test the decoy localization candidates next to the possible sites."""

    def test_adjacent_decoy_candidates(self):
        """Decoys on the priority residue (anywhere) and on adjacent residues next to possible sites."""
        search = pd.DataFrame(
            {
                "RAW_FILE": ["a", "a"],
                "SCAN_NUMBER": [1, 1],
                "SEQUENCE": ["KGWPSAK"] * 2,
                "PRECURSOR_CHARGE": [2, 2],
                "MODIFIED_SEQUENCE": ["KGWPS[UNIMOD:21]AK", "KGWPS[UNIMOD:21]AK"],
            }
        )
        result = add_adjacent_decoy_candidates(search, {21: ["S", "T", "Y"]}, ["G"], ["A", "G", "W"], max_distance=2)
        decoys = sorted(result["MODIFIED_SEQUENCE"].iloc[2:])
        # G: priority residue (position 2), A: next to S (position 6), W: 2 residues from S (position 3)
        self.assertEqual(decoys, sorted(["KG[UNIMOD:21]WPSAK", "KGWPSA[UNIMOD:21]K", "KGW[UNIMOD:21]PSAK"]))
        self.assertEqual(len(result), 5)

    def test_adjacent_decoys_keep_other_sites(self):
        """With two phosphorylations, one moves to the decoy site and the other stays on a possible site."""
        search = pd.DataFrame(
            {
                "RAW_FILE": ["a"],
                "SCAN_NUMBER": [1],
                "SEQUENCE": ["PSATK"],
                "PRECURSOR_CHARGE": [2],
                "MODIFIED_SEQUENCE": ["PS[UNIMOD:21]AT[UNIMOD:21]K"],
            }
        )
        result = add_adjacent_decoy_candidates(search, {21: ["S", "T"]}, [], ["A"], max_distance=1)
        self.assertEqual(
            sorted(result["MODIFIED_SEQUENCE"].iloc[1:]),
            sorted(["PS[UNIMOD:21]A[UNIMOD:21]TK", "PSA[UNIMOD:21]T[UNIMOD:21]K"]),
        )


class TestLocalizationGroups(unittest.TestCase):
    """Test the groups used to estimate the false localization rate."""

    def test_groups_and_neighbors(self):
        """Residue of the called site and whether a possible site is directly next to it."""
        peptides = pd.Series(["K.PES[UNIMOD:21]TK.R", "-.PES[UNIMOD:21]PTK.-", "-.AY[UNIMOD:21]PK.-"])
        self.assertEqual(localization_groups(peptides, r"\[UNIMOD:21\]").tolist(), ["S", "S", "Y"])
        self.assertEqual(
            neighbor_classes(peptides, r"\[UNIMOD:21\]", "STY").tolist(), ["neighbor", "no neighbor", "no neighbor"]
        )


class TestDecoyWeights(unittest.TestCase):
    """Test the expected false localizations of a decoy win."""

    def test_count_weights(self):
        """Test a decoy win stands for (real candidates - 1) / decoy candidates false localizations."""
        from oktoberfest.rescore.localization import _decoy_weights

        psms = pd.DataFrame(
            {
                "filename": ["f"] * 5 + ["g"] * 2,
                "ScanNr": [1] * 5 + [2] * 2,
                "peptide": [
                    "_.S[UNIMOD:21]TAK._",
                    "_.ST[UNIMOD:21]AK._",
                    "_.STA[UNIMOD:21]K._",
                    "_.STAA[UNIMOD:21]._",
                    "_.STA[UNIMOD:21]A._",
                    "_.S[UNIMOD:21]GK._",
                    "_.SG[UNIMOD:21]K._",
                ],
                "score_prob": [1.0, 0.5, 0.2, 0.1, 0.3, 1.0, 0.4],
            }
        )
        decoys = psms["peptide"].str.contains(r"[AG]\[UNIMOD:21\]")
        weights = _decoy_weights(psms, decoys, r"\[UNIMOD:21\]")
        # spectrum 1: 2 real, 3 decoy candidates; spectrum 2: 1 real, 1 decoy candidate
        np.testing.assert_allclose(weights.loc[("f", 1)], 1 / 3)
        np.testing.assert_allclose(weights.loc[("g", 2)], 0.0)


class TestHybridFlr(unittest.TestCase):
    """Test the combination of the decoy and the model based FLR."""

    def test_hybrid_localization_flr(self):
        """Test the hybrid FLR is max(model high, min(decoy, model low)), with model FLRs per group."""
        from oktoberfest.rescore import hybrid_localization_flr

        psms = pd.DataFrame(
            {
                "q-value": [0.001, 0.001, 0.001, 0.001, 0.5],
                "loc_prob": [0.9, 0.8, 0.7, 0.9, 0.9],
                "localization_flr": [0.0, 0.0, 0.05, 0.0, np.nan],
                "model_error_high": [0.02, 0.04, 0.3, 0.1, 0.0],
                "model_error_low": [0.01, 0.02, 0.2, 0.05, 0.0],
                "model_group": ["S neighbor", "S neighbor", "S neighbor", "Y no neighbor", "S neighbor"],
            }
        )
        result = hybrid_localization_flr(psms.copy(), "model_group")
        # group "S neighbor" (the PSM above the PSM FDR is not used): running means of the sorted errors
        high = [0.02, 0.03, 0.12]
        low = [0.01, 0.015, 0.23 / 3]
        decoy = [0.0, 0.0, 0.05]
        expected = [max(h, min(d, lo)) for h, d, lo in zip(high, decoy, low, strict=True)]
        np.testing.assert_allclose(result["localization_flr"].iloc[:3], expected)
        # its own group: model FLRs 0.1 / 0.05, decoy 0 -> max(0.1, min(0, 0.05))
        np.testing.assert_allclose(result["localization_flr"].iloc[3], 0.1)
        self.assertTrue(np.isnan(result["localization_flr"].iloc[4]))
        np.testing.assert_allclose(result["localization_flr_decoy"].iloc[:4], [0.0, 0.0, 0.05, 0.0])


class TestScanNumbers(unittest.TestCase):
    """Test the scan numbers parsed from percolator PSM ids."""

    def test_numeric_raw_file_names(self):
        """Test raw files named by numbers only (read as int by pandas) give the right scan numbers."""
        from oktoberfest.rescore.localization import _scan_numbers

        psms = pd.DataFrame({"PSMId": ["12-345-PEPS[UNIMOD:21]K-2-1", "file-7-PEPK-2-1"], "filename": [12, "file"]})
        self.assertEqual(_scan_numbers(psms), [345, 7])


class TestAmbiguityCheck(unittest.TestCase):
    """Test the site-determining ions, the ambiguity flags and the ranking of ambiguous localizations (V1)."""

    def test_differing_ions(self):
        """Test only the b / y ions with a different set of PTMs differ, up to fragment charge min(3, charge)."""
        from oktoberfest.rescore.localization import _differing_ions

        # sites 4 vs 5 of PEPSTK: b4 and y2 (= complement of b4) differ, at fragment charges 1 and 2
        differs = _differing_ions("PEPS[UNIMOD:21]TK", "PEPST[UNIMOD:21]K", 2).reshape(5, 6)
        expected = np.zeros((5, 6), dtype=bool)
        expected[1, :2] = True  # y2+1, y2+2
        expected[3, 3:5] = True  # b4+1, b4+2
        np.testing.assert_array_equal(differs, expected)
        # same localization: nothing differs; other modifications are no site-determining ions
        self.assertFalse(
            _differing_ions("[UNIMOD:1]-M[UNIMOD:35]S[UNIMOD:21]K", "[UNIMOD:1]-M[UNIMOD:35]S[UNIMOD:21]K", 3).any()
        )

    def test_ambiguity_flags(self):
        """Test calls without an observed site-determining ion, or whose observed ones favour B, are ambiguous."""
        from unittest import mock

        from oktoberfest.rescore import localization as loc

        a, b = "PEPS[UNIMOD:21]TK", "PEPST[UNIMOD:21]K"
        localization = pd.DataFrame(
            {
                "PSMId": [f"f-{scan}-{seq}-2-1" for scan in (1, 2, 3, 4) for seq in (a, b)],
                "filename": "f",
                "ScanNr": [1, 1, 2, 2, 3, 3, 4, 4],
                "peptide": [f"_.{seq}._" for _ in range(4) for seq in (a, b)],
                "score": [2.0, 1.0] * 4,
                "loc_prob": [0.7, 0.3] * 4,
            }
        )
        # one row per (scan, sequence): possible but unobserved (1e-7) everywhere, observed b4+1 / y2+1 (index 3 * 6 + 3
        # of A and 1 * 6 of B) where given
        raw = np.full((8, 30), 1e-7)
        raw[2, 21] = 0.5  # scan 2: b4+1 of A observed -> evidence for A
        raw[4, 21], raw[5, 6] = 0.2, 0.6  # scan 3: B's y2+1 stronger than A's b4+1 -> favours B
        rows = {(scan, seq, 2): 2 * (scan - 1) + j for scan in (1, 2, 3, 4) for j, seq in enumerate((a, b))}
        del rows[(4, b, 2)]  # scan 4: no stored annotation of B -> not flagged
        with mock.patch.object(loc, "_stored_annotation", return_value=(rows, np.zeros_like(raw), raw)):
            result = loc.ambiguity_flags(localization, Path("."))
        calls = result[result["loc_prob"] == 0.7].set_index("ScanNr")
        self.assertEqual(calls["localization_ambiguous"].tolist(), [True, False, True, False])
        self.assertEqual(calls.loc[1, "localization_runner_up"], f"_.{b}._")
        self.assertFalse(result.loc[result["loc_prob"] == 0.3, "localization_ambiguous"].any())

    def _psms(self):
        return pd.DataFrame(
            {
                "q-value": [0.001] * 4,
                "loc_prob": [0.9, 0.8, 0.7, 0.6],
                "loc_prob_delta": [0.8, 0.6, 0.4, 0.2],
                "decoy_weight": [1.0, 0.0, 0.0, 0.0],
                "model_error_high": [0.01, 0.02, 0.03, 0.04],
                "model_error_low": [0.005, 0.01, 0.015, 0.02],
                "model_group": "S neighbor",
            }
        )

    def test_ambiguous_ranked_last(self):
        """Test ambiguous calls are ranked after the others by both estimates, with their own errors (V1)."""
        from oktoberfest.rescore import estimate_localization_flr, hybrid_localization_flr

        psms = self._psms().assign(localization_ambiguous=[True, False, False, False], decoy_weight=0.0)
        result = hybrid_localization_flr(estimate_localization_flr(psms), "model_group")
        np.testing.assert_allclose(result["localization_flr_decoy"], 0.0)
        # model: order 2, 3, 4, 1 -> running means of the original errors, as q-values; decoy 0 -> the high estimate
        np.testing.assert_allclose(result["localization_flr"], [0.025, 0.02, 0.025, 0.025])

    def test_ambiguous_decoy_floor(self):
        """Test the FLR of ambiguous calls is at least the decoy FLR of the ambiguous calls alone, pooled (E1a)."""
        from oktoberfest.rescore import estimate_localization_flr, hybrid_localization_flr

        psms = pd.concat([self._psms(), self._psms().assign(model_group="Y no neighbor", decoy_weight=0.0)])
        psms = psms.reset_index(drop=True).assign(localization_ambiguous=[True, False, False, True] * 2)
        result = hybrid_localization_flr(estimate_localization_flr(psms), "model_group")
        # decoy over all: order 2, 3, 2', 3', then the ambiguous 1, 1', 4, 4' -> one decoy win in 8, as q-values
        np.testing.assert_allclose(result["localization_flr_decoy"], [0.125, 0.0, 0.0, 0.125] * 2)
        # model per group: order 2, 3, 1, 4 -> running means 0.02, 0.025, 0.02, 0.025 of the errors, as q-values
        high = [0.02, 0.02, 0.02, 0.025]
        low = [0.01, 0.01, 0.01, 0.0125]
        v1 = [max(h, min(d, lo)) for h, d, lo in zip(high * 2, result["localization_flr_decoy"], low * 2, strict=True)]
        # ambiguous calls alone: 1, 1', 4, 4' by loc_prob_delta -> cumulative weights 1, 1, 1, 1 / rank -> q 1 / 4
        floor = [0.25, 0.0, 0.0, 0.25] * 2
        np.testing.assert_allclose(result["localization_flr"], np.maximum(v1, floor))

    def test_without_ambiguous_calls_unchanged(self):
        """Test no ambiguous call gives bit-identical estimates to the estimates without the ambiguity check."""
        from oktoberfest.rescore import estimate_localization_flr, hybrid_localization_flr

        off = hybrid_localization_flr(estimate_localization_flr(self._psms()), "model_group")
        on = hybrid_localization_flr(
            estimate_localization_flr(self._psms().assign(localization_ambiguous=False)), "model_group"
        )
        for column in ("localization_flr", "localization_flr_decoy"):
            np.testing.assert_array_equal(off[column].to_numpy(), on[column].to_numpy())


class TestLocalizationScore(unittest.TestCase):
    """Test the localization score calibrated on the estimated FLR."""

    def test_anchors(self):
        """Test 1% FLR maps to 0.75, 5% to 0.5, linear in between and to the ends, NaN without an FLR."""
        from oktoberfest.rescore import calibrate_localization_score

        psms = pd.DataFrame({"localization_flr": [0.0, 0.005, 0.01, 0.03, 0.05, 0.525, 1.0, np.nan]})
        score = calibrate_localization_score(psms)["localization_score"]
        np.testing.assert_allclose(score, [1.0, 0.875, 0.75, 0.625, 0.5, 0.25, 0.0, np.nan])

    def test_score_cut_equals_flr_cut(self):
        """Test a score cut selects exactly the localizations at the matching estimated FLR, also unordered."""
        from oktoberfest.rescore import calibrate_localization_score

        rng = np.random.default_rng(1)
        psms = pd.DataFrame({"localization_flr": rng.uniform(0, 0.2, 1000), "loc_prob_delta": rng.uniform(0, 1, 1000)})
        score = calibrate_localization_score(psms)["localization_score"]
        np.testing.assert_array_equal(score >= 0.75, psms["localization_flr"] <= 0.01)
        np.testing.assert_array_equal(score >= 0.5, psms["localization_flr"] <= 0.05)
