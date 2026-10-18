import unittest

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
