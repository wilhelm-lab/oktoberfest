import unittest
from pathlib import Path

import pandas as pd
from numpy.testing import assert_almost_equal

from oktoberfest.data import Spectra
from oktoberfest.data.spectra import FragmentType
from oktoberfest.pr import Predictor


class TestTMTProsit(unittest.TestCase):
    """Test class for TMT model predictions."""

    def test_prosit_tmt(self):
        """Test retrieval of predictions from prosit tmt models via koina."""
        meta_df = pd.read_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        var = Spectra._gen_vars_df()
        library = Spectra(obs=meta_df, var=var)
        library.strings_to_categoricals()

        intensity_predictor = Predictor.from_koina(
            model_name="Prosit_2020_intensity_TMT",
            server_url="koina.wilhelmlab.org:443",
            ssl=True,
            targets=["intensities", "annotation"]
        )
        intensity_predictor.predict_intensities(data=library)

        irt_predictor = Predictor.from_koina(
            model_name="Prosit_2020_irt_TMT",
            server_url="koina.wilhelmlab.org:443",
            ssl=True,
        )
        irt_predictor.predict_rt(data=library)

        library_expected = Spectra.from_hdf5(Path(__file__).parent / "data" / "predictions" / "library_output.h5ad.gz")

        assert_almost_equal(
            library.get_matrix(FragmentType.PRED)[0].toarray(),
            library_expected.get_matrix(FragmentType.PRED)[0].toarray(),
            decimal=6,
        )
        # explicitly set this to int64 as the default on windows is int32
        # which causes the type check to fail on windows
        library.var["num"] = library.var["num"].astype("int64")
        library.var["charge"] = library.var["charge"].astype("int64")

        pd.testing.assert_frame_equal(library.obs, library_expected.obs)
        pd.testing.assert_frame_equal(library.var, library_expected.var)

