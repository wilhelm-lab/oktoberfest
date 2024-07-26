import contextlib
import logging
import multiprocessing as mp
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, Union

import numpy as np
from spectrum_fundamentals.constants import ALPHABET, ANNOTATION
from spectrum_fundamentals.mod_string import parse_modstrings
from spectrum_io.file.parquet import write_file

sys.path.append("/cmnfs/home/students/j.schlensok/dlomix/bmpc_shared_scripts/oktoberfest_interface")
# Suppress unnecessary feature processor info from DLomix
with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):  # noqa: E402
    from oktoberfest_interface import (
        download_model_from_github,
        load_keras_model,
        process_dataset,
    )

from ..data.spectra import Spectra  # noqa: E402

logger = logging.getLogger(__name__)

ANNOTATIONS = [f"{ion_type}{pos}+{charge}".encode() for ion_type, charge, pos in list(zip(*ANNOTATION))]


class DLomix:
    """A class for interacting with DLomix models locally for inference."""

    def __init__(self, model_type: str, model_path: Union[Path, str], output_path: Path):
        """Initialize a DLomix predictor from name or path of pre-loaded weights.

        Weights for a baseline intensity predictor can also be downloaded by specifying model_path=\"baseline\"

       :param model_type: Type of model (intensity or irt)
       :param model_path: Path of pre-trained PrositIntensityPredictor, or \"baseline\"
       :param output_path: Directory to save processed data for predictor to for reuse

       :raises NotImplementedError: if a retention time predictor is requested
       :raises ValueError: if a name other than \"baseline\" is provided as model_path instead of a path
        """
        self.model_type = model_type
        self.model_path = model_path
        self.output_path = output_path

        if model_type == "irt":
            logger.exception(NotImplementedError())

        self.output_name = "intensities"

        if model_path == "baseline":
            with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
                model_path = download_model_from_github()
        elif isinstance(model_path, str):
            raise ValueError("Non-path model names other than 'baseline' are not supported for local prediction")

        logger.info(f"Loading model weights from {model_path}")
        self.model = load_keras_model(str(model_path))

    def predict(self, data: Spectra, dataset_name: str, temporary: bool = False) -> Dict[str, np.ndarray]:
        """Create predictions for dataset using Keras model.

        :param data: spectral library to predict features for
        :param dataset_name: Name of the dataset for storing processed files for DLomix
        :param temporary: Whether to keep the processed files on the disk after they've been used for inference

        :return: a dictionary containing predicted features (key: feature type) and a mask of the ion annotations of
            the predicted feature matrix (key: 'annotation')

        """
        # TODO create wrapper for folder with Parquet file + ion type & modification files
        logger.info(f"Pre-processing dataset {dataset_name}")

        data_dir = self.output_path / dataset_name

        parquet_path = data_dir / "processed_dataset.parquet"
        ion_type_file = data_dir / "ion_types.txt"
        modification_file = data_dir / "modifications.txt"

        create_new_parquet = True

        if parquet_path.exists():
            if ion_type_file.exists() and modification_file.exists():
                logger.info(f"Parquet file {parquet_path} already exists, reusing it")
                create_new_parquet = False
                ion_types = np.loadtxt(ion_type_file, dtype=object).tolist()
                modifications = np.loadtxt(modification_file, dtype=object).tolist()
            else:
                logger.info(
                    f"""Existing Parquet file found at {parquet_path} but missing ion type/modification info, so it
                     will be overwritten"""
                )

        data_dir.mkdir(exist_ok=True, parents=True)

        if create_new_parquet:
            processed_data = data.assemble_df_for_parquet(include_intensities=False)
            write_file(processed_data, parquet_path)
            ion_types = data.uns["ion_types"].tolist()
            modifications = sorted(
                list(
                    {
                        token
                        for peptide in parse_modstrings(processed_data["modified_sequence"].tolist(), ALPHABET)
                        for token in peptide
                    }
                ),
                key=lambda token: ALPHABET[token],
            )

        if not temporary:
            ion_type_file.write_text("\n".join(ion_types))
            modification_file.write_text("\n".join(modifications))

        with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds = process_dataset(
                str(parquet_path),
                self.model,
                ion_types=ion_types,
                modifications=modifications,
                label_column=None,
                test_ratio=1,
                val_ratio=0,
            )

        # TODO set batch size based on available memory
        # TODO improve Keras progress bar:
        # - solve glitching with multi-threading
        # - possibly write custom callback that only returns progress information so we can display it ourselves
        preds = self.model.predict(ds.tensor_inference_data)

        if temporary:
            parquet_path.unlink()
            data_dir.rmdir()

        return {self.output_name: preds, "annotation": np.tile(np.array(ANNOTATIONS), (preds.shape[0], 1))}

    @staticmethod
    def initialize_tensorflow(gpu_number: int = 0) -> None:
        """
        Set some enviroment variables before importing TensorFlow.

        - only use one GPU (otherwise TensorFlow will reserve all available ones)
        - use all available cores
        - change logging level

        :param gpu_number: Number of GPU to use
        """
        logging.getLogger("tensorflow").setLevel(logging.WARNING)
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_number)

        n_threads = mp.cpu_count()

        os.environ["OP_NUM_THREADS"] = str(n_threads)
        os.environ["TF_NUM_INTRAOP_THREADS"] = str(n_threads)
        os.environ["TF_NUM_INTEROP_THREADS"] = str(n_threads)
