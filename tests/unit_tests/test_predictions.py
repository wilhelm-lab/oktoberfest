import unittest
from pathlib import Path

import pandas as pd
from numpy.testing import assert_almost_equal

from oktoberfest.data import Spectra
from oktoberfest.data.spectra import FragmentType
from oktoberfest.pr import predict


class TestTMTProsit(unittest.TestCase):
    """Test class for TMT model predictions."""

    def test_prosit_tmt(self):
        """Test retrieval of predictions from prosit tmt models via koina."""
        meta_df = pd.read_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        var = Spectra._gen_vars_df()
        library = Spectra(obs=meta_df, var=var)
        library.strings_to_categoricals()
        pred_intensities = predict(
            library.obs,
            model_name="Prosit_2020_intensity_TMT",
            server_url="koina.wilhelmlab.org:443",
            ssl=True,
            targets=["intensities", "annotation"],
        )
        pred_irt = predict(
            library.obs, model_name="Prosit_2020_irt_TMT", server_url="koina.wilhelmlab.org:443", ssl=True
        )

        library.add_matrix(pred_intensities["intensities"], FragmentType.PRED)
        library.add_column(pred_irt["irt"].squeeze(), name="PREDICTED_IRT")

        library_expected = Spectra.from_hdf5(Path(__file__).parent / "data" / "predictions" / "library_output.h5ad.gz")

        assert_almost_equal(
            library.get_matrix(FragmentType.PRED)[0].toarray(),
            library_expected.get_matrix(FragmentType.PRED)[0].toarray(),
            decimal=6,
            verbose=True,
        )
        pd.testing.assert_frame_equal(library.obs, library_expected.obs)
        pd.testing.assert_frame_equal(library.var, library_expected.var)
