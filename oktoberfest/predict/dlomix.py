"""DLOmix local predictor using InferencePipeline."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..data import Spectra

from dlomix.reports.postprocessing import normalize_intensity_predictions

logger = logging.getLogger(__name__)


class DLOmixLocal:
    """Local DLOmix predictor using InferencePipeline.

    Uses a pre-trained DLOmix model bundled with its preprocessor.
    Similar interface to Koina for intensity prediction in rescoring.
    """

    def __init__(self, pipeline):
        """
        Initialize DLOmixLocal predictor.

        :param pipeline: Loaded InferencePipeline instance
        """
        self.pipeline = pipeline

    def predict(
        self,
        data: Any,
        **kwargs,
    ) -> dict[str, np.ndarray]:
        """
        Predict fragment ion intensities using DLOmix InferencePipeline.

        :param data: Spectra object with observations (modified_sequence, precursor_charge, collision_energy)
        :param kwargs: Additional arguments (unused, for API compatibility)
        :return: Dictionary with predictions:
                 - 'intensities': predicted fragment intensities (n_psms, n_fragments)
                 - 'annotation': fragment ion annotations
                 - 'mz': m/z values (zeros)
        """
        n_input_psms = data.n_obs if hasattr(data, 'n_obs') else data.shape[0]

        logger.info(f"Running DLOmix predictions for {n_input_psms} PSMs")

        # Build input dataframe from Spectra.obs with required columns
        input_df = self._build_input_df(data)

        # Run predictions
        predictions = self.pipeline.predict(input_df)

        # Get expected fragment annotations from Spectra var_names
        fragment_annotations = np.array(data.var_names) if hasattr(data, 'var_names') else None

        # Process output with normalization
        result = self._process_output(predictions, data=data, input_df=input_df, fragment_annotations=fragment_annotations)
        logger.info("Successfully predicted intensities")

        return result

    def _build_input_df(self, data: Any) -> pd.DataFrame:
        """
        Build input dataframe for InferencePipeline from Spectra object.

        Extracts required model features and transforms collision_energy by normalizing by 100.
        Preprocessor handles sequence encoding and one-hot encoding of precursor_charge.

        :param data: Spectra object
        :return: DataFrame with required columns
        """
        seq_col = self.pipeline.preprocessor.sequence_column
        model_features = self.pipeline.preprocessor.model_features

        col_map = {c.lower(): c for c in data.obs.columns}
        input_data = {}

        # Get sequence column
        seq_col_lower = seq_col.lower()
        if seq_col_lower in col_map:
            input_data[seq_col] = data.obs[col_map[seq_col_lower]].values
        else:
            raise ValueError(f"Sequence column '{seq_col}' not found in Spectra.obs")

        # Transform model features
        for feat in model_features:
            # Map feature name to column name (handle collision_energy_aligned_normed -> collision_energy)
            if feat == "collision_energy_aligned_normed":
                col_name = col_map.get("collision_energy")
                if col_name is None:
                    raise ValueError("collision_energy column not found in Spectra.obs")
                values = data.obs[col_name].values.astype(float)
                input_data[feat] = (values / 100.0).astype(np.float32)
            else:
                feat_lower = feat.lower()
                if feat_lower not in col_map:
                    raise ValueError(f"Required feature '{feat}' not found in Spectra.obs")
                col_name = col_map[feat_lower]
                input_data[feat] = data.obs[col_name].values

        return pd.DataFrame(input_data)

    def _process_output(
        self,
        intensities: np.ndarray,
        data: Any = None,
        input_df: pd.DataFrame = None,
        fragment_annotations: Optional[np.ndarray] = None,
    ) -> dict[str, np.ndarray]:
        """
        Process model output into Oktoberfest format with normalization.

        :param intensities: Model predictions (numpy array, shape n_psms x n_fragments)
        :param data: Spectra object with sequence information
        :param input_df: Input DataFrame with precursor_charge
        :param fragment_annotations: Fragment ion annotations
        :return: Dictionary with 'intensities', 'annotation', 'mz'
        """
        # Ensure 2D: (n_psms, n_fragments)
        if len(intensities.shape) == 1:
            intensities = intensities.reshape(-1, 1)

        n_fragments = intensities.shape[1]
        n_psms = intensities.shape[0]

        # Apply DLOmix normalization
        if data is not None and input_df is not None:
            col_map = {c.lower(): c for c in data.obs.columns}
            seq_col = col_map.get("modified_sequence") or col_map.get("sequence")

            charges = input_df["precursor_charge"].values[:n_psms].astype(int)
            onehot = np.zeros((len(charges), 6), dtype=np.float32)
            for i, c in enumerate(charges):
                if 1 <= c <= 6:
                    onehot[i, c - 1] = 1.0

            norm_df = pd.DataFrame({
                "sequences": data.obs[seq_col].values[:n_psms],
                "intensities_pred": [list(row) for row in intensities],
                "precursor_charge_onehot": [list(row) for row in onehot],
            })

            norm_df = normalize_intensity_predictions(norm_df, compute_spectral_angle=False)
            intensities = np.array(norm_df["intensities_pred"].tolist(), dtype=np.float32)

        # Use provided fragment annotations or create placeholder indices
        if fragment_annotations is not None and len(fragment_annotations) == n_fragments:
            annotation = np.array([fragment_annotations])
        else:
            annotation = np.array([np.arange(1, n_fragments + 1, dtype=object).astype(str)])

        # m/z values: set to 1.0 (DLOmix doesn't predict m/z)
        mz = np.ones((n_psms, n_fragments), dtype=np.float32)

        return {
            'intensities': intensities.astype(np.float32),
            'annotation': annotation,
            'mz': mz,
        }
