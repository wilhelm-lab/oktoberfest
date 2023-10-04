import unittest
from pathlib import Path

import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.pr import grpc_predict


class TestTMTProsit(unittest.TestCase):
    """Test class for TMT model predictions."""

    def test_prosit_tmt(self):
        """Test retrieval of predictions from prosit tmt models via koina."""
        library = Spectra.from_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        grpc_predict(
            library=library,
            url="koina.proteomicsdb.org:443",
            intensity_model="Prosit_2020_intensity_TMT",
            irt_model="Prosit_2020_irt_TMT",
            ssl=True,
            alignment=False,
            job_type="",
        )

        expected_df = pd.read_csv(Path(__file__).parent / "data" / "predictions" / "library_output.csv")
        sparse_cols = [col for col in library.spectra_data.columns if col.startswith("INTENSITY_PRED")]
        for sparse_col in sparse_cols:
            expected_df[sparse_col] = expected_df[sparse_col].astype(library.spectra_data[sparse_col].dtype)
        expected_df["PREDICTED_IRT"] = expected_df["PREDICTED_IRT"].astype(library.spectra_data["PREDICTED_IRT"].dtype)

        pd.testing.assert_frame_equal(library.spectra_data, expected_df)
