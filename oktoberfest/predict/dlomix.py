"""DLOmix local predictor for fine-tuned models integration with Oktoberfest."""

from __future__ import annotations

import json
import logging
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd
import tensorflow as tf

if TYPE_CHECKING:
    from ..data import Spectra

logger = logging.getLogger(__name__)


def extract_alphabet_from_model(model_path: str | Path) -> Optional[dict]:
    """
    Extract the alphabet dictionary from a Keras model file.

    Keras .keras files are ZIP archives containing a config.json with model configuration.
    This function extracts the alphabet dictionary that was used during training.

    :param model_path: Path to the .keras model file
    :return: Dictionary mapping amino acids to integer codes, or None if not found
    """
    try:
        model_path = Path(model_path)
        if not model_path.exists():
            logger.warning(f"Model file not found: {model_path}")
            return None

        with zipfile.ZipFile(str(model_path), 'r') as zip_file:
            if 'config.json' not in zip_file.namelist():
                logger.warning("No config.json found in model file")
                return None

            with zip_file.open('config.json') as f:
                config = json.load(f)
                alphabet = config.get('config', {}).get('alphabet')

                if alphabet:
                    logger.info(f"Extracted alphabet from model: {len(alphabet)} tokens")
                    return alphabet
                else:
                    logger.warning("No alphabet found in model config")
                    return None
    except Exception as e:
        logger.warning(f"Failed to extract alphabet from model: {e}")
        return None


class DLOmixLocal:
    """Local DLOmix predictor for fine-tuned models.

    Similar to Koina, but runs predictions locally using a fine-tuned DLOmix model.
    Accepts input as Spectra, pd.DataFrame, or tensorflow dataset.
    Automatically extracts alphabet from model config for proper sequence encoding.
    """

    def __init__(self, model_path: str | Path, parquet_path: Optional[str | Path] = None, extract_alphabet: bool = True):
        """
        Initialize DLOmixLocal predictor.

        :param model_path: Path to the .keras model weights file
        :param parquet_path: Path to parquet file with pre-processed features
        :param extract_alphabet: If True, automatically extract alphabet from model config
        :raises FileNotFoundError: if model file does not exist
        :raises ImportError: if tensorflow is not installed
        """
        self.model_path = Path(model_path)
        if not self.model_path.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")

        self.parquet_path = Path(parquet_path) if parquet_path else None
        self.parquet_data = None
        self.mz_data = None
        if self.parquet_path and self.parquet_path.exists():
            logger.info(f"Loading pre-processed features from {parquet_path}")
            self.parquet_data = pd.read_parquet(str(self.parquet_path))
            # Extract mz_raw if available
            if 'mz_raw' in self.parquet_data.columns:
                self.mz_data = self.parquet_data['mz_raw'].values

        self.alphabet = None
        if extract_alphabet:
            self.alphabet = extract_alphabet_from_model(self.model_path)
            if self.alphabet:
                logger.info(f"Extracted alphabet from model ({len(self.alphabet)} tokens)")

        self.model = self._load_model()
        logger.info(f"Loaded DLOmix model from {model_path}")

    def _load_model(self) -> Any:
        """Load TensorFlow Keras model."""
        return tf.keras.models.load_model(str(self.model_path))

    def predict(
        self,
        data: Any,
        **kwargs,
    ) -> dict[str, np.ndarray]:
        """
        Predict fragment ion intensities using fine-tuned DLOmix model.

        Uses dlomix dataset processing for proper sequence encoding.

        :param data: Input data for prediction (Spectra object with fragment annotations)
        :param kwargs: Additional arguments (unused, for API compatibility)
        :return: Dictionary with predictions:
                 - 'intensities': predicted fragment intensities (n_psms, n_fragments)
                 - 'annotation': fragment ion annotations
                 - 'mz': fragment m/z values
        """
        logger.info("Running DLOmix predictions")

        if self.parquet_path is None:
            raise ValueError("DLOmixLocal requires parquet_path for prediction")

        # Extract fragment annotations from Spectra object if available
        fragment_annotations = None
        if hasattr(data, 'var_names'):
            fragment_annotations = np.array(data.var_names)

        # Get actual PSM count from Spectra object
        n_input_psms = None
        if hasattr(data, 'n_obs'):
            n_input_psms = data.n_obs
        elif hasattr(data, 'shape'):
            n_input_psms = data.shape[0]

        # Use dlomix dataset to handle sequence encoding properly
        try:
            from dlomix.data import FragmentIonIntensityDataset
        except ImportError:
            raise ImportError("dlomix is required for sequence encoding")

        parquet_source = str(self.parquet_path)

        # Filter parquet to match batch if full file is larger
        parquet_to_use = parquet_source
        mz_data_trimmed = self.mz_data

        if self.parquet_data is not None and len(self.parquet_data) > n_input_psms:
            # Try to match by scan number if available
            if hasattr(data, 'obs') and 'SCAN' in data.obs.columns:
                batch_scans = set(data.obs['SCAN'].values)
                batch_parquet = self.parquet_data[self.parquet_data['scan_number'].isin(batch_scans)].copy()
                logger.info(f"Filtering parquet by scan number: {len(self.parquet_data)} → {len(batch_parquet)} PSMs")
            else:
                # Fallback: just take first N rows
                batch_parquet = self.parquet_data.iloc[:n_input_psms].copy()
                logger.info(f"Filtering parquet by index: {len(self.parquet_data)} → {len(batch_parquet)} PSMs")

            # Write temp parquet file
            import tempfile
            temp_dir = tempfile.gettempdir()
            parquet_to_use = f"{temp_dir}/batch_data_{id(batch_parquet)}.parquet"
            batch_parquet.to_parquet(parquet_to_use)
            mz_data_trimmed = batch_parquet['mz_raw'].values if 'mz_raw' in batch_parquet.columns else None
        else:
            logger.info(f"Loading all {len(self.parquet_data) if self.parquet_data is not None else '?'} PSMs from parquet")

        dataset = FragmentIonIntensityDataset(
            test_data_source=parquet_to_use,
            data_format="parquet",
            sequence_column="modified_sequence",
            label_column="intensities_raw",
            encoding_scheme="unmod",
            max_seq_len=30,
            model_features=[
                "precursor_charge_onehot",
                "collision_energy_aligned_normed",
            ],
            alphabet=self.alphabet,
            dataset_type="tf",
        )

        # Get test data
        if hasattr(dataset, 'get_tf_dataset'):
            test_data = dataset.get_tf_dataset('test')
        elif hasattr(dataset, 'tensor_test_data'):
            test_data = dataset.tensor_test_data
        else:
            test_data = dataset.hf_dataset

        predictions = self.model.predict(test_data, verbose=0)

        result = self._process_output(predictions, mz_data=mz_data_trimmed, fragment_annotations=fragment_annotations)
        logger.info("Successfully predicted intensities")

        return result

    def _process_output(
        self,
        intensities: Any,
        mz_data: Optional[np.ndarray] = None,
        fragment_annotations: Optional[np.ndarray] = None,
    ) -> dict[str, np.ndarray]:
        """
        Process model output into Oktoberfest format.

        :param intensities: Model predictions (TensorFlow tensor or numpy array)
        :param mz_data: Optional m/z values from parquet file
        :param fragment_annotations: Fragment ion annotations (e.g., "b1+1", "y1+1")
        :return: Dictionary with 'intensities', 'annotation', 'mz'
        """
        # Convert to numpy if needed
        if hasattr(intensities, 'numpy'):
            intensities = intensities.numpy()
        else:
            intensities = np.array(intensities)

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

        # Use mz_data if available, otherwise zeros
        if mz_data is not None:
            # Handle array of arrays case
            if len(mz_data) > 0 and isinstance(mz_data[0], np.ndarray):
                mz_array = np.stack(mz_data, dtype=np.float32)
                # Trim or pad to match n_psms
                if mz_array.shape[0] < n_psms:
                    mz = np.zeros((n_psms, n_fragments), dtype=np.float32)
                    mz[:mz_array.shape[0]] = mz_array
                else:
                    mz = mz_array[:n_psms]
            else:
                mz_data = np.array(mz_data, dtype=np.float32)
                # If mz_data is already 2D, trim or pad to match n_psms
                if mz_data.ndim == 2 and mz_data.shape[1] == n_fragments:
                    if mz_data.shape[0] < n_psms:
                        mz = np.zeros((n_psms, n_fragments), dtype=np.float32)
                        mz[:mz_data.shape[0]] = mz_data
                    else:
                        mz = mz_data[:n_psms]
                # Otherwise create per-fragment mz (broadcast to all PSMs)
                elif mz_data.ndim == 1:
                    mz = np.tile(mz_data, (n_psms, 1))
                else:
                    mz = mz_data[:n_psms]
        else:
            mz = np.zeros((n_psms, n_fragments), dtype=np.float32)

        return {
            'intensities': intensities,
            'annotation': annotation,
            'mz': mz,
        }
