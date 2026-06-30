"""DLOmix local predictor using InferencePipeline."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..data import Spectra

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

        :param data: Spectra object with observations (modified_sequence, precursor_charge_onehot, collision_energy_aligned_normed)
        :param kwargs: Additional arguments (unused, for API compatibility)
        :return: Dictionary with predictions:
                 - 'intensities': predicted fragment intensities (n_psms, n_fragments)
                 - 'annotation': fragment ion annotations
                 - 'mz': m/z values (zeros)
        """
        # Get actual PSM count and fragment annotations
        n_input_psms = data.n_obs if hasattr(data, 'n_obs') else data.shape[0]
        fragment_annotations = np.array(data.var_names) if hasattr(data, 'var_names') else None

        logger.info(f"Running DLOmix predictions for {n_input_psms} PSMs")

        # Build input dataframe from Spectra.obs with required columns
        input_df = self._build_input_df(data)

        # Run predictions
        predictions = self.pipeline.predict(input_df)

        # Process output
        result = self._process_output(predictions, fragment_annotations=fragment_annotations)
        logger.info("Successfully predicted intensities")

        return result

    def _build_input_df(self, data: Any) -> pd.DataFrame:
        """
        Build input dataframe for InferencePipeline from Spectra object.

        Extracts and transforms: modified_sequence, precursor_charge_onehot, collision_energy_aligned_normed.

        :param data: Spectra object
        :return: DataFrame with required columns
        """
        seq_col = self.pipeline.preprocessor.sequence_column
        model_features = self.pipeline.preprocessor.model_features

        # Create lowercase column mapping for case-insensitive lookup
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
            if feat == "precursor_charge_onehot":
                # Convert PRECURSOR_CHARGE (1-6) to one-hot encoding
                charge_col = col_map.get("precursor_charge") or col_map.get("charge")
                if charge_col is None:
                    raise ValueError("PRECURSOR_CHARGE column not found")
                charges = data.obs[charge_col].values.astype(int)
                # One-hot encode: charges 1-6 → 6-dimensional vector
                onehot = np.zeros((len(charges), 6), dtype=np.float32)
                for i, c in enumerate(charges):
                    if 1 <= c <= 6:
                        onehot[i, c - 1] = 1.0
                input_data[feat] = onehot.tolist()
            elif feat == "collision_energy_aligned_normed":
                # Normalize collision energy by dividing by 100
                ce_col = col_map.get("collision_energy")
                if ce_col is None:
                    raise ValueError("COLLISION_ENERGY column not found")
                ce_values = data.obs[ce_col].values.astype(float)
                ce_normalized = ce_values / 100.0
                input_data[feat] = ce_normalized.astype(np.float32)
            else:
                col_lower = feat.lower()
                if col_lower in col_map:
                    input_data[feat] = data.obs[col_map[col_lower]].values
                else:
                    raise ValueError(f"Required feature '{feat}' not found in Spectra.obs")

        return pd.DataFrame(input_data)

    def _process_output(
        self,
        intensities: np.ndarray,
        fragment_annotations: Optional[np.ndarray] = None,
    ) -> dict[str, np.ndarray]:
        """
        Process model output into Oktoberfest format.

        :param intensities: Model predictions (numpy array, shape n_psms x n_fragments)
        :param fragment_annotations: Fragment ion annotations
        :return: Dictionary with 'intensities', 'annotation', 'mz'
        """
        # Ensure 2D: (n_psms, n_fragments)
        if len(intensities.shape) == 1:
            intensities = intensities.reshape(-1, 1)

        n_fragments = intensities.shape[1]
        n_psms = intensities.shape[0]

        # Use provided fragment annotations or create placeholder indices
        if fragment_annotations is not None and len(fragment_annotations) == n_fragments:
            annotation = np.array([fragment_annotations])
        else:
            annotation = np.array([np.arange(1, n_fragments + 1, dtype=object).astype(str)])

        # m/z values default to zeros (DLOmix doesn't predict m/z)
        mz = np.zeros((n_psms, n_fragments), dtype=np.float32)

        return {
            'intensities': intensities.astype(np.float32),
            'annotation': annotation,
            'mz': mz,
        }
