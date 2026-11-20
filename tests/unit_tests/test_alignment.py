import unittest

import numpy as np
import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.data.spectra import FragmentType
from oktoberfest.predict.alignment import _alignment, _prepare_alignment_df

CE_RANGE = (30, 33)


def _make_library(rows: list[dict]) -> Spectra:
    obs = pd.DataFrame(rows)
    obs.index = obs.index.astype(str)
    return Spectra(obs=obs, var=Spectra._gen_vars_df())


def _psm(scan: int, score: float, ce: float = 30.0, charge: int = 2, fragmentation: str = "HCD", reverse=False):
    return {
        "RAW_FILE": "run_a",
        "SCAN_NUMBER": scan,
        "SCORE": score,
        "COLLISION_ENERGY": ce,
        "PRECURSOR_CHARGE": charge,
        "FRAGMENTATION": fragmentation,
        "REVERSE": reverse,
    }


class TestPrepareAlignmentDf(unittest.TestCase):
    """Test the selection of PSMs used for collision energy calibration."""

    def test_selects_hcd_targets_and_best_psm_per_scan(self):
        """Test decoys, non-HCD spectra and lower scoring PSMs of the same scan are excluded."""
        library = _make_library(
            [
                _psm(scan=1, score=10),
                _psm(scan=1, score=50),  # better PSM of scan 1
                _psm(scan=2, score=40, reverse=True),
                _psm(scan=3, score=30, fragmentation="CID"),
                _psm(scan=4, score=20),
            ]
        )

        alignment_library = _prepare_alignment_df(library, ce_range=CE_RANGE)

        n_ces = CE_RANGE[1] - CE_RANGE[0]
        self.assertEqual(alignment_library.n_obs, 2 * n_ces)
        selected = alignment_library.obs.drop_duplicates("SCAN_NUMBER")
        self.assertEqual(selected["SCAN_NUMBER"].tolist(), [1, 4])
        self.assertEqual(selected["SCORE"].tolist(), [50, 20])
        self.assertTrue((alignment_library.obs["FRAGMENTATION"] == "HCD").all())
        self.assertFalse(alignment_library.obs["REVERSE"].any())

    def test_propagates_each_psm_to_every_ce(self):
        """Test every selected PSM is repeated for each CE and keeps its original CE."""
        library = _make_library([_psm(scan=1, score=50, ce=27.0), _psm(scan=2, score=40, ce=35.0)])

        alignment_library = _prepare_alignment_df(library, ce_range=CE_RANGE)

        obs = alignment_library.obs
        self.assertEqual(len(obs), 6)
        for scan, orig_ce in [(1, 27.0), (2, 35.0)]:
            psm_rows = obs[obs["SCAN_NUMBER"] == scan]
            self.assertEqual(sorted(psm_rows["COLLISION_ENERGY"].tolist()), list(range(*CE_RANGE)))
            self.assertTrue((psm_rows["ORIG_COLLISION_ENERGY"] == orig_ce).all())
        pd.testing.assert_frame_equal(alignment_library.var, library.var)

    def test_top_n_per_ce_or_per_charge(self):
        """Test the top PSMs are selected per original CE for DDA and per precursor charge for DIA."""
        # 25 PSMs at CE 30 and 25 at CE 35, all charge 2; the CE 35 PSMs score higher
        library = _make_library(
            [_psm(scan=i, score=i, ce=30.0) for i in range(25)]
            + [_psm(scan=i, score=i, ce=35.0) for i in range(25, 50)]
        )

        # crosslinked peptides use the top 20 PSMs per group
        per_ce = _prepare_alignment_df(library, ce_range=(30, 31), group_by_charge=False, xl=True)
        per_charge = _prepare_alignment_df(library, ce_range=(30, 31), group_by_charge=True, xl=True)

        # top 20 of each original CE
        self.assertEqual(sorted(per_ce.obs["SCAN_NUMBER"].tolist()), list(range(5, 25)) + list(range(30, 50)))
        # top 20 of the single charge group, independent of the original CE
        self.assertEqual(sorted(per_charge.obs["SCAN_NUMBER"].tolist()), list(range(30, 50)))

    def test_top_n_limits(self):
        """Test at most 1000 PSMs per group are used for linear and 20 for crosslinked peptides."""
        library = _make_library([_psm(scan=i, score=i) for i in range(1, 1003)])

        linear = _prepare_alignment_df(library, ce_range=(30, 31))
        xl = _prepare_alignment_df(library, ce_range=(30, 31), xl=True)

        self.assertEqual(linear.n_obs, 1000)
        self.assertEqual(xl.n_obs, 20)
        # the highest scoring PSMs are kept
        self.assertEqual(xl.obs["SCAN_NUMBER"].min(), 983)


class TestAlignment(unittest.TestCase):
    """Test the spectral angle calculation used for collision energy calibration."""

    def setUp(self):  # noqa: D102
        self.library = _make_library([_psm(scan=1, score=1), _psm(scan=2, score=1)])
        self.n_vars = self.library.n_vars

    def test_spectral_angle(self):
        """Test identical spectra get a spectral angle of 1 and spectra without shared peaks of 0."""
        raw = np.full((2, self.n_vars), -1.0)
        pred = np.full((2, self.n_vars), -1.0)
        raw[:, :3] = [1.0, 0.5, 0.0]
        pred[0, :3] = [1.0, 0.5, 0.0]  # identical
        pred[1, :3] = [0.0, 0.0, 1.0]  # orthogonal
        self.library.add_intensities_without_mapping(raw, FragmentType.RAW)
        self.library.add_intensities_without_mapping(pred, FragmentType.PRED)

        _alignment(self.library)

        np.testing.assert_allclose(self.library.obs["SPECTRAL_ANGLE"], [1.0, 0.0], atol=1e-6)

    def test_spectral_angle_xl(self):
        """Test the crosslinked spectral angle is the mean of the angles of both peptides."""
        identical = np.full((2, self.n_vars), -1.0)
        identical[:, :2] = [1.0, 0.5]
        orthogonal = np.full((2, self.n_vars), -1.0)
        orthogonal[:, :2] = [0.0, 1.0]
        raw_b = np.full((2, self.n_vars), -1.0)
        raw_b[:, :2] = [1.0, 0.0]
        self.library.add_intensities_without_mapping(identical.copy(), FragmentType.RAW_A)
        self.library.add_intensities_without_mapping(identical.copy(), FragmentType.PRED_A)
        self.library.add_intensities_without_mapping(raw_b, FragmentType.RAW_B)
        self.library.add_intensities_without_mapping(orthogonal, FragmentType.PRED_B)

        _alignment(self.library, xl=True)

        np.testing.assert_allclose(self.library.obs["SPECTRAL_ANGLE_A"], [1.0, 1.0], atol=1e-6)
        np.testing.assert_allclose(self.library.obs["SPECTRAL_ANGLE_B"], [0.0, 0.0], atol=1e-6)
        np.testing.assert_allclose(self.library.obs["SPECTRAL_ANGLE"], [0.5, 0.5], atol=1e-6)
