import logging
import multiprocessing as mp
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import spectrum_fundamentals.constants as c
import spectrum_io.file.parquet as parquet
from spectrum_fundamentals.fragments import format_fragment_ion_annotation, generate_fragment_ion_annotations
from spectrum_fundamentals.mod_string import parse_modstrings

from oktoberfest.data import Spectra
from oktoberfest.utils import Config
from oktoberfest.utils.logging import mute_stdout

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

with mute_stdout():
    from dlomix.interface import download_model_from_github, load_keras_model, process_dataset, save_keras_model
    from dlomix.refinement_transfer_learning.automatic_rl_tl import AutomaticRlTlTraining, AutomaticRlTlTrainingConfig

gpus = tf.config.list_physical_devices("GPU")
for device in gpus:
    tf.config.experimental.set_memory_growth(device, True)

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

logger = logging.getLogger(__name__)
logging.getLogger("dlomix").setLevel(logging.INFO)

ANNOTATIONS = [f"{ion_type}{pos}+{charge}".encode() for ion_type, charge, pos in list(zip(*c.ANNOTATION))]
OPTIMAL_ION_TYPE_ORDER = [
    "y",
    "b",
    "x",
    "z",
    "z_r",
    "a",
    "c",
]  # y > b > rest so that intensity predictor can re-use weights


def _download_baseline_model(model_path: Path) -> None:
    logger.info(f"Downloading baseline intensity predictor to {model_path}")
    downloaded_model_path = Path(download_model_from_github())
    downloaded_model_path.rename(model_path)


def refine_intensity_predictor(
    baseline_model_path: Path,
    libraries: List[Spectra],
    config: Config,
    data_directory: Path,
    result_directory: Path,
    dataset_name: str,
    model_name: str,
    download_new_baseline_model: bool = False,
) -> None:
    """Perform refinement/transfer learning on a baseline intensity predictor.

    :param baseline_model_path: Path of baseline model to refine
    :param libraries: Spectral libraries to use as training data for refinement
    :param config: Config containing refinement learning options
    :param data_directory: Directory to save processed dataset and refined model to. The Parquet and ion_type and
        modification metadata files for the processed dataset will be stored in `<output_directory>/<dataset_name>/`,
        and the refined model as `<output_directory>/<model_name>.keras`
    :param result_directory: Directory to save CSV logs & report notebook to
    :param dataset_name: Name of dataset
    :param model_name: Name of refined model
    :param download_new_baseline_model: Whether to download a new baseline model from GitHub to the specified path
    """
    gpus = tf.config.list_physical_devices("GPU")
    if len(gpus) == 0:
        logger.warning(
            """TensorFlow could not detect GPU devices to use for refinement learning, so it will run on CPU and take
            much longer"""
        )

    if config.include_original_sequences:
        for library in libraries:
            library.obs["MODIFIED_SEQUENCE_RAW"] = library.obs["MODIFIED_SEQUENCE"]
            additional_columns = ["SEQUENCE", "MODIFIED_SEQUENCE_RAW"]
    else:
        additional_columns = []

    parquet_path, ion_types, modifications = create_dlomix_dataset(
        libraries,
        data_directory / dataset_name,
        include_additional_columns=additional_columns,
        andromeda_score_threshold=config.andromeda_score_threshold,
        num_duplicates=config.num_duplicates,
    )

    if download_new_baseline_model:
        _download_baseline_model(baseline_model_path)
    baseline_model = load_keras_model(str(baseline_model_path))

    model_path = data_directory / (model_name + ".keras")
    if model_path.exists():
        logger.info(f"Found existing refined model at {model_path}, re-using it")
        return

    additional_columns = [column.lower() for column in additional_columns]
    logger.info("Pre-processing dataset for refinement learning")
    ds = process_dataset(
        str(parquet_path),
        baseline_model,
        ion_types=ion_types,
        modifications=modifications,
        label_column="intensities_raw",
        val_ratio=0.2,
        batch_size=config.training_batch_size,
        additional_columns=additional_columns,
    )

    training_config = AutomaticRlTlTrainingConfig(
        dataset=ds,
        baseline_model=baseline_model,
        improve_further=config.improve_further,
        use_wandb=config.use_wandb,
        wandb_project=config.wandb_project,
        wandb_tags=config.wandb_tags,
        results_log=str(result_directory),
    )

    trainer = AutomaticRlTlTraining(training_config)
    refined_model = trainer.train()

    save_keras_model(refined_model, str(model_path))


def create_dlomix_dataset(
    libraries: List[Spectra],
    output_dir: Path,
    include_additional_columns: Optional[List[str]] = None,
    andromeda_score_threshold: Optional[float] = None,
    num_duplicates: Optional[int] = None,
) -> Tuple[Path, List[str], List[str]]:
    """Transform one or multiple spectra into Parquet file that can be used by DLomix.

    Processes spectral libraries into DLomix-compatible format and detects fragment ion types and peptide modifications
        present in the dataset, then writes the dataset to `output_dir` as `processed_dataset.parquet` and the lists of
        ion types and modifications as `ion_types.txt` and `modifications.txt`.

    :param libraries: Spectral libraries to include
    :param output_dir: Directory to save processed dataset to
    :param include_additional_columns: additional columns to keep in the dataset
    :param andromeda_score_threshold: Andromeda score cutoff for peptides included in output
    :param num_duplicates: Number of (sequence, charge, collision energy) duplicates to keep in output

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
    modifications = {
        modification
        for spectra in libraries
        for peptide in parse_modstrings(spectra.obs["MODIFIED_SEQUENCE"].tolist(), alphabet)
        for modification in peptide
    }
    token = max(alphabet.values()) + 1
    for new_modification in modifications - set(alphabet):
        alphabet[new_modification] = token
        token += 1
    modifications = sorted(list(modifications), key=lambda modification: alphabet[modification])
    logger.debug(f"Detected modifications in dataset: {modifications}")
    modification_file.write_text("\n".join(modifications))

    ion_types = sorted(
        list({ion_type for spectra in libraries for ion_type in spectra.uns["ion_types"]}),
        key=lambda ion: OPTIMAL_ION_TYPE_ORDER.index(ion),
    )
    logger.debug(f"Detected ion types in dataset: {ion_types}")
    ion_type_file.write_text("\n".join(ion_types))

    all_data = Spectra(ad.concat([spectra for spectra in libraries]))
    processed_data = all_data.preprocess_for_machine_learning(
        ion_type_order=ion_types,
        include_additional_columns=include_additional_columns,
        andromeda_score_threshold=andromeda_score_threshold,
        num_duplicates=num_duplicates,
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
        annotations = np.tile(np.array(fragment_ion_order, dtype=object), (preds.shape[0], 1))

        if not keep_dataset:
            shutil.rmtree(self.output_path / dataset_name)

        return {self.output_name: preds, "annotation": annotations}
