import unittest
from pathlib import Path

import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.data.spectra import FragmentType
from oktoberfest.pr import predict


class TestTMTProsit(unittest.TestCase):
    """Test class for TMT model predictions."""

    def test_prosit_tmt(self):
        """Test retrieval of predictions from prosit tmt models via koina."""
        library = Spectra.from_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        input_data = library.spectra_data

        pred_intensities = predict(
            input_data,
            model_name="Prosit_2020_intensity_TMT",
            server_url="koina.proteomicsdb.org:443",
            ssl=True,
            targets=["intensities", "annotation"],
        )
        pred_irt = predict(
            input_data, model_name="Prosit_2020_irt_TMT", server_url="koina.proteomicsdb.org:443", ssl=True
        )

        library.add_matrix(pd.Series(pred_intensities["intensities"].tolist(), name="intensities"), FragmentType.PRED)
        library.add_column(pred_irt["irt"], name="PREDICTED_IRT")

        expected_df = pd.read_csv(Path(__file__).parent / "data" / "predictions" / "library_output.csv")
        sparse_cols = library.get_matrix(FragmentType.PRED)[1]
        for sparse_col in range(0, len(sparse_cols)):
            expected_df[sparse_cols[sparse_col]] = expected_df[sparse_cols[sparse_col]].astype(
                library.spectra_data.layers["pred_int"][:, sparse_col].dtype
            )
        expected_df["PREDICTED_IRT"] = expected_df["PREDICTED_IRT"].astype(
            library.spectra_data.obs["PREDICTED_IRT"].dtype
        )

        # pd.testing.assert_frame_equal(library.spectra_daa, expected_df)

    def test_failing_koina(self):
        """Test koina with input data that does not fit to the model to trigger exception handling."""
        library = Spectra.from_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        input_data = library.spectra_data

        self.assertRaises(
            Exception,
            predict,
            input_data,
            model_name="Prosit_2020_intensity_HCD",
            server_url="koina.proteomicsdb.org:443",
            ssl=True,
            targets=["intensities", "annotation"],
        )
