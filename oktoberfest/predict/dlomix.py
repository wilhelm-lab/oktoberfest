import contextlib
import logging
import multiprocessing as mp
import os
import shutil
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
import spectrum_io.file.parquet as parquet
from spectrum_fundamentals.fragments import generate_fragment_ion_annotations
from spectrum_fundamentals.mod_string import parse_modstrings

from oktoberfest.data import Spectra

# TODO maybe move this to DLomix::interface::__init__.py?
# Set some enviroment variables before importing TensorFlow.
# 1. change logging level
logging.getLogger("tensorflow").setLevel(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

# 2. use all available cores
n_threads = mp.cpu_count()
os.environ["OP_NUM_THREADS"] = str(n_threads)
os.environ["TF_NUM_INTRAOP_THREADS"] = str(n_threads)
os.environ["TF_NUM_INTEROP_THREADS"] = str(n_threads)

import tensorflow as tf
from dlomix.interface import download_model_from_github, load_keras_model, process_dataset, save_keras_model
from dlomix.refinement_transfer_learning.automatic_rl_tl import AutomaticRlTlTraining, AutomaticRlTlTrainingConfig

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

logger = logging.getLogger(__name__)

ANNOTATIONS = [f"{ion_type}{pos}+{charge}".encode() for ion_type, charge, pos in list(zip(*c.ANNOTATION))]
OPTIMAL_ION_TYPE_ORDER = ["y", "b", "x", "z", "a", "c"]  # y > b > rest so that intensity predictor can re-use weights


@contextlib.contextmanager
def mute_stdout(ignore_warnings: bool = False):
    """Mute print statements and user warnings from DLomix."""
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
        if ignore_warnings:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                yield
        else:
            yield


def _download_baseline_model(model_path: Path) -> None:
    logger.info(f"Downloading baseline intensity predictor to {model_path}")
    downloaded_model_path = Path(download_model_from_github())
    downloaded_model_path.rename(model_path)


def refine_intensity_predictor(
    baseline_model_path: Path,
    spectra: List[Spectra],
    data_directory: Path,
    result_directory: Path,
    dataset_name: str,
    model_name: str,
    download_new_baseline_model: bool = False,
    batch_size: int = 1024,
    additional_columns: Optional[List[str]] = None,
    available_gpus: Optional[List[int]] = None,
    use_wandb: bool = False,
    wandb_project: Optional[str] = None,
    wandb_tags: Optional[List[str]] = None,
) -> None:
    """Perform refinement/transfer learning on a baseline intensity predictor.

    :param baseline_model_path: Path of baseline model to refine
    :param spectra: Spectral libraries to use as training data for refinement
    :param data_directory: Directory to save processed dataset and refined model to. The Parquet and ion_type and
        modification metadata files for the processed dataset will be stored in `<output_directory>/<dataset_name>/`,
        and the refined model as `<output_directory>/<model_name>.keras`
    :param result_directory: Directory to save CSV logs & report notebook to.
    :param dataset_name: Name of dataset
    :param model_name: Name of refined model
    :param download_new_baseline_model: Whether to download a new baseline model from GitHub to the specified path.
    :param batch_size: Batch size to use for training
    :param additional_columns: Additional columns to keep in DLomix dataset for downstream analyis
    :param available_gpus: Indices of GPUs to use for training
    :param use_wandb: Whether to use WandB to log training
    :param wandb_project: Name of WandB project to save run to
    :param wandb_tags: Tags to assing to WandB run
    """
    gpus = tf.config.list_physical_devices("GPU")
    if len(gpus) == 0:
        logger.warning(
            """TensorFlow could not detect GPU devices to use for refinement learning, so it will run on CPU and take
            much longer"""
        )
    else:
        for device in gpus:
            tf.config.experimental.set_memory_growth(device, True)
        if available_gpus:
            tf.config.set_visible_devices([gpus[i] for i in available_gpus], "GPU")

    parquet_path, ion_types, modifications = create_dlomix_dataset(
        spectra, data_directory / dataset_name, include_additional_columns=additional_columns
    )

    if download_new_baseline_model:
        _download_baseline_model(baseline_model_path)
    baseline_model = load_keras_model(str(baseline_model_path))

    model_path = data_directory / (model_name + ".keras")
    if model_path.exists():
        logger.info(f"Found existing refined model at {model_path}, re-using it")
        return

    if additional_columns:
        additional_columns = [column_name.lower() for column_name in additional_columns]
    else:
        additional_columns = []

    logger.info("Pre-processing dataset for refinement learning")
    with mute_stdout(ignore_warnings=True):
        ds = process_dataset(
            str(parquet_path),
            baseline_model,
            ion_types=ion_types,
            modifications=modifications,
            label_column="intensities_raw",
            val_ratio=0.2,
            batch_size=batch_size,
            additional_columns=additional_columns,
        )

    config = AutomaticRlTlTrainingConfig(
        dataset=ds,
        baseline_model=baseline_model,
        use_wandb=use_wandb,
        wandb_project=wandb_project,
        wandb_tags=wandb_tags,
        results_log=str(result_directory),
    )

    trainer = AutomaticRlTlTraining(config)
    refined_model = trainer.train()

    save_keras_model(refined_model, str(model_path))


def create_dlomix_dataset(
    spectra: List[Spectra], output_dir: Path, include_additional_columns: Optional[List[str]] = None
) -> Tuple[Path, List[str], List[str]]:
    """Transform one or multiple spectra into Parquet file that can be used by DLomix.

    Processes spectral libraries into DLomix-compatible format and detects fragment ion types and peptide modifications
        present in the dataset, then writes the dataset to `output_dir` as `processed_dataset.parquet` and the lists of
        ion types and modifications as `ion_types.txt` and `modifications.txt`.

    :param spectra: Spectral libraries to include
    :param output_dir: Directory to save processed dataset to
    :param include_additional_columns: additional columns to keep in the dataset

    :returns:
        - path of saved Parquet file
        - a list of ion types in it
        - a list of modifications in it (in the form of modstring tokens from spectrum_fundamentals)
    """
    # TODO create wrapper for folder with Parquet file + ion type & modification files

    output_dir.mkdir(exist_ok=True, parents=True)
    parquet_path = output_dir / "processed_dataset.parquet"
    ion_type_file = output_dir / "ion_types.txt"
    modification_file = output_dir / "modifications.txt"

    # Re-use existing dataset if all necessary files are there
    if parquet_path.exists() and ion_type_file.exists() and modification_file.exists():
        logger.info(f"Re-using saved dataset from {output_dir}")
        ion_types = np.loadtxt(ion_type_file, dtype=object).tolist()
        modifications = np.loadtxt(modification_file, dtype=object).tolist()
        return parquet_path, ion_types, modifications

    # Otherwise regenerate dataset because ion type metadata isn't stored in the Parquet file
    if not include_additional_columns:
        include_additional_columns = []
    include_additional_columns += c.SHARED_DATA_COLUMNS

    alphabet = c.ALPHABET
    val = [max(alphabet.values()) + 1]

    # TODO fix custom modification parsing in fundamentals
    def get_alphabet_value(token, alphabet, val):
        if token not in alphabet:
            alphabet[token] = val[0]
            val[0] += 1
        return alphabet[token]

    modifications = sorted(
        list(
            {
                token
                for spectrum in spectra
                for peptide in parse_modstrings(spectrum.obs["MODIFIED_SEQUENCE"].tolist(), alphabet)
                for token in peptide
            }
        ),
        key=lambda token: get_alphabet_value(token, alphabet, val),
    )
    logger.debug(f"Detected modifications in dataset: {modifications}")
    modification_file.write_text("\n".join(modifications))

    ion_types = sorted(
        list({ion_type for spectrum in spectra for ion_type in spectrum.uns["ion_types"]}),
        key=lambda ion: OPTIMAL_ION_TYPE_ORDER.index(ion),
    )
    logger.debug(f"Detected ion types in dataset: {ion_types}")
    ion_type_file.write_text("\n".join(ion_types))

    processed_data = pd.concat(
        [
            spectrum.preprocess_for_machine_learning(
                ion_type_order=ion_types, include_additional_columns=include_additional_columns
            )
            for spectrum in spectra
        ]
    )
    if not parquet_path.exists():
        logger.info(f"Saving DLomix dataset to {parquet_path}")
        parquet.write_file(processed_data, parquet_path)

    return parquet_path, ion_types, modifications


class DLomix:
    """A class for interacting with DLomix models locally for inference."""

    def __init__(self, model_type: str, model_path: Path, output_path: Path, batch_size: int, download: bool = False):
        """Initialize a DLomix predictor from name or path of pre-loaded weights.

        :param model_type: Type of model (intensity or irt)
        :param model_path: Path of model file
        :param output_path: Directory to save processed data for predictor to for reuse
        :param batch_size: Batch size to use for inference
        :param download: Whether to download a new baseline model from GitHub to the model_path

        :raises NotImplementedError: if a retention time predictor is requested
        """
        self.model_type = model_type
        self.output_path = output_path
        self.batch_size = batch_size

        if model_type == "irt":
            raise NotImplementedError("Local prediction not implemented for iRT prediction")

        self.output_name = "intensities"
        if download:
            _download_baseline_model(model_path)
        self.model = load_keras_model(str(model_path))

    def predict(self, data: Spectra, dataset_name: str, keep_dataset: bool = True) -> Dict[str, np.ndarray]:
        """Create predictions for dataset using Keras model.

        :param data: spectral library to predict features for
        :param dataset_name: Name of the dataset for storing processed files for DLomix
        :param keep_dataset: Whether to keep or discard the pre-processed dataset after inference

        :return: a dictionary containing predicted features (key: feature type) and a mask of the ion annotations of
            the predicted feature matrix (key: 'annotation')

        """
        # TODO reuse training dataset if doing transfer learning (load subset based on RAW_FILE column)
        parquet_path, ion_types, modifications = create_dlomix_dataset([data], self.output_path / dataset_name)
        with mute_stdout(ignore_warnings=True):
            ds = process_dataset(
                str(parquet_path),
                self.model,
                ion_types=ion_types,
                modifications=modifications,
                label_column=None,
                test_ratio=1,
                val_ratio=0,
                batch_size=self.batch_size,
            )

        # TODO improve Keras progress bar:
        # - solve glitching with multi-threading
        # - possibly write custom callback that only returns progress information so we can display it ourselves
        preds = self.model.predict(ds.tensor_inference_data)
        fragment_ion_order = [
            format_fragment_ion_annotation(ann)
            for ann in generate_fragment_ion_annotations(ion_types, order=("position", "ion_type", "charge"))
        ]

        if not keep_dataset:
            shutil.rmtree(self.output_path / dataset_name)

        return {self.output_name: preds, "annotation": fragment_ion_order}
