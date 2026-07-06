import unittest

import pandas as pd

from .. import DATA_DIR
from oktoberfest import rescore as re
from oktoberfest.data import Spectra


class TestRescoring(unittest.TestCase):
    """Test class for rescoring functions."""

    def test_generate_features(self):
        """Test feature generation with internal results and additional columns."""
        # load test library containing 100 psms
        library = Spectra.from_hdf5(DATA_DIR / "library" / "library100.hdf5")

        # all additional columns
        re.generate_features(
            library=library,
            search_type="rescore",
            output_file=DATA_DIR / "library" / "rescore_all.tab",
            additional_columns="all",
            all_features=False,
            regression_method="spline",
        )
        expected_all = pd.read_csv(DATA_DIR / "library" / "expected_rescore_all.tab", sep="\t")
        created_all = pd.read_csv(DATA_DIR / "library" / "rescore_all.tab", sep="\t")
        pd.testing.assert_frame_equal(expected_all, created_all)

        # no additional columns
        re.generate_features(
            library=library,
            search_type="rescore",
            output_file=DATA_DIR / "library" / "rescore_none.tab",
            additional_columns="none",
            all_features=False,
            regression_method="spline",
        )
        expected_none = pd.read_csv(DATA_DIR / "library" / "expected_rescore_none.tab", sep="\t")
        created_none = pd.read_csv(DATA_DIR / "library" / "rescore_none.tab", sep="\t")
        pd.testing.assert_frame_equal(expected_none, created_none)

        # list of additional columns
        re.generate_features(
            library=library,
            search_type="rescore",
            output_file=DATA_DIR / "library" / "rescore_list.tab",
            additional_columns=["A"],
            all_features=False,
            regression_method="spline",
        )
        expected_list = pd.read_csv(DATA_DIR / "library" / "expected_rescore_list.tab", sep="\t")
        created_list = pd.read_csv(DATA_DIR / "library" / "rescore_list.tab", sep="\t")
        pd.testing.assert_frame_equal(expected_list, created_list)

        # load test library containing 150 psms for multifrag
        library = Spectra.from_hdf5(DATA_DIR / "library" / "library150_multifrag.hdf5")

        # all additional columns
        re.generate_features(
            library=library,
            search_type="rescore",
            output_file=DATA_DIR / "library" / "rescore_all_multifrag.tab",
            additional_columns=None,
            all_features=False,
            task="multifrag",
            regression_method="spline",
            featured_ions=["C", "z"],  # c,z ions
        )

        expected_all = pd.read_csv(DATA_DIR / "library" / "expected_rescore_multifrag.tab", sep="\t")
        created_all = pd.read_csv(DATA_DIR / "library" / "rescore_all_multifrag.tab", sep="\t")
        pd.testing.assert_frame_equal(expected_all, created_all)
