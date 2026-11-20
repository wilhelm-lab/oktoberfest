import importlib.util
import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.predict.dlomix import DLOmixLocal

DLOMIX_INSTALLED = importlib.util.find_spec("dlomix") is not None


def _make_pipeline(model_features: list[str], sequence_column: str = "modified_sequence") -> MagicMock:
    pipeline = MagicMock()
    pipeline.preprocessor = SimpleNamespace(sequence_column=sequence_column, model_features=model_features)
    return pipeline


def _make_library(fragmentation: list[str], charges: list[int], sequences: list[str] = None) -> Spectra:
    n = len(charges)
    obs = pd.DataFrame(
        {
            "MODIFIED_SEQUENCE": sequences or ["PEPTIDEK"] * n,
            "PRECURSOR_CHARGE": charges,
            "COLLISION_ENERGY": [30.0] * n,
            "FRAGMENTATION": fragmentation,
        }
    )
    var = Spectra._gen_vars_df()
    library = Spectra(obs=obs, var=var)
    library.var_names = var.index
    return library


class TestDLOmixInput(unittest.TestCase):
    """Test building the InferencePipeline input from Spectra."""

    def test_build_input_df_derived_features(self):
        """Test CE normalization and one-hot encoding of charge and fragmentation method."""
        predictor = DLOmixLocal(
            _make_pipeline(["collision_energy_aligned_normed", "precursor_charge_onehot", "method_onehot"])
        )
        library = _make_library(fragmentation=["HCD", "CID"], charges=[2, 3])

        input_df = predictor._build_input_df(library)

        self.assertEqual(input_df["modified_sequence"].tolist(), ["PEPTIDEK", "PEPTIDEK"])
        np.testing.assert_allclose(input_df["collision_energy_aligned_normed"], [0.3, 0.3])
        self.assertEqual(input_df["collision_energy_aligned_normed"].dtype, np.float32)
        np.testing.assert_array_equal(
            np.stack(input_df["precursor_charge_onehot"]), [[0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]]
        )
        # one-hot over (CID, HCD)
        np.testing.assert_array_equal(np.stack(input_df["method_onehot"]), [[0, 1], [1, 0]])

    def test_build_input_df_passes_through_obs_columns(self):
        """Test model features available in obs are passed through (case-insensitive)."""
        predictor = DLOmixLocal(_make_pipeline(["precursor_charge"]))
        library = _make_library(fragmentation=["HCD"], charges=[4])

        input_df = predictor._build_input_df(library)

        self.assertEqual(input_df["precursor_charge"].tolist(), [4])

    def test_build_input_df_rejects_unsupported_fragmentation(self):
        """Test method_onehot fails for fragmentation methods the model was not trained on."""
        predictor = DLOmixLocal(_make_pipeline(["method_onehot"]))
        for fragmentation in ["ETD", "unknown"]:
            with self.subTest(fragmentation=fragmentation):
                library = _make_library(fragmentation=["HCD", fragmentation], charges=[2, 2])
                with self.assertRaisesRegex(ValueError, "CID and HCD"):
                    predictor._build_input_df(library)

    def test_build_input_df_missing_columns(self):
        """Test clear errors when the sequence column or a model feature is missing."""
        library = _make_library(fragmentation=["HCD"], charges=[2])
        with self.assertRaisesRegex(ValueError, "Sequence column"):
            DLOmixLocal(_make_pipeline([], sequence_column="peptide_sequence"))._build_input_df(library)
        with self.assertRaisesRegex(ValueError, "instrument_type"):
            DLOmixLocal(_make_pipeline(["instrument_type"]))._build_input_df(library)


@unittest.skipUnless(DLOMIX_INSTALLED, "requires the dlomix extra")
class TestDLOmixPredict(unittest.TestCase):
    """Test predictions and post-processing of the DLOmix output."""

    def test_predict_normalizes_and_masks_impossible_ions(self):
        """Test predictions are normalized to the base peak and impossible ions are set to -1."""
        pipeline = _make_pipeline(["collision_energy_aligned_normed", "precursor_charge_onehot"])
        rng = np.random.default_rng(42)
        raw_predictions = rng.uniform(-0.1, 2.0, size=(2, 174)).astype(np.float32)
        pipeline.predict.return_value = raw_predictions
        library = _make_library(fragmentation=["HCD", "HCD"], charges=[1, 3], sequences=["PEPTIDEK", "PEPTIDEKAAR"])

        result = DLOmixLocal(pipeline).predict(library)

        pipeline.predict.assert_called_once()
        intensities = result["intensities"]
        self.assertEqual(intensities.shape, (2, 174))
        self.assertEqual(intensities.dtype, np.float32)
        np.testing.assert_allclose(intensities.max(axis=1), [1.0, 1.0], rtol=1e-6)
        valid = intensities[intensities != -1]
        self.assertTrue(((valid >= 0) & (valid <= 1)).all())

        annotations = np.array(library.var_names)
        fragment_num = library.var["num"].to_numpy()
        fragment_charge = library.var["charge"].to_numpy()
        # fragments longer than the peptide cannot exist
        np.testing.assert_array_equal(intensities[0, fragment_num >= 8], -1)
        np.testing.assert_array_equal(intensities[1, fragment_num >= 11], -1)
        # fragment charges above the precursor charge cannot exist for precursor charge 1
        np.testing.assert_array_equal(intensities[0, fragment_charge > 1], -1)
        self.assertTrue((intensities[1, (fragment_charge == 3) & (fragment_num < 11)] != -1).all())

        np.testing.assert_array_equal(result["annotation"], [annotations])
        np.testing.assert_array_equal(result["mz"], np.ones((2, 174)))
