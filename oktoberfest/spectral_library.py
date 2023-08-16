import logging
import re
from math import ceil
from multiprocessing import current_process
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from spectrum_io.file import csv
from spectrum_io.spectral_library import digest
from tqdm.auto import tqdm
from tritonclient.grpc import InferenceServerClient, InferInput, InferRequestedOutput

from .data.spectra import FragmentType, Spectra
from .utils.config import Config

logger = logging.getLogger(__name__)


def infer_predictions(
    triton_client: InferenceServerClient,
    model: str,
    input_data: Dict[str, Tuple[np.ndarray, str]],
    outputs: List[str],
    batch_size: int,
) -> Dict[str, np.ndarray]:
    """
    Infer predictions from a triton client.

    :param triton_client: An inference client using grpc
    :param model: a model that is recognized by the server specified in the triton client
    :param input_data: a dictionary that contains the input names (key) for the specific model
        and a tuple of the input_data as a numpy array of shape [:, 1] and the dtype recognized
        by the triton client (value).
    :param outputs: a list of output names for the specific model
    :param batch_size: the number of elements from the input_data that should be provided to the
        triton client at once
    :return: a dictionary containing the predictions (values) for the given outputs (keys)
    """
    num_spec = len(input_data[list(input_data)[0]][0])
    predictions: Dict[str, List[np.ndarray]] = {output: [] for output in outputs}

    n_batches = ceil(num_spec / batch_size)
    process_identity = current_process()._identity
    if len(process_identity) > 0:
        position = process_identity[0]
    else:
        position = 0

    with tqdm(
        total=n_batches,
        position=position,
        desc=f"Inferring predictions for {num_spec} spectra with batch site {batch_size}",
        leave=True,
    ) as progress:
        for i in range(0, n_batches):
            progress.update(1)
            # logger.info(f"Predicting batch {i+1}/{n_batches}.")
            infer_inputs = []
            for input_key, (data, dtype) in input_data.items():
                batch_data = data[i * batch_size : (i + 1) * batch_size]
                infer_input = InferInput(input_key, batch_data.shape, dtype)
                infer_input.set_data_from_numpy(batch_data)
                infer_inputs.append(infer_input)

            infer_outputs = [InferRequestedOutput(output) for output in outputs]

            prediction = triton_client.infer(model, inputs=infer_inputs, outputs=infer_outputs)

            for output in outputs:
                predictions[output].append(prediction.as_numpy(output))

    return {key: np.vstack(value) for key, value in predictions.items()}


def parse_fragment_labels(spectra_labels: np.ndarray, precursor_charges: np.ndarray, seq_lengths: np.ndarray):
    """Uses regex to parse labels."""
    pattern = rb"([y|b])([0-9]{1,2})\+([1-3])"
    fragment_types = []
    fragment_numbers = []
    fragment_charges = []
    for spectrum_labels in spectra_labels:
        types = []
        numbers = []
        charges = []
        for label in spectrum_labels:
            match = re.match(pattern, label)
            if match:
                groups = match.groups()
                types.append(groups[0].decode())
                numbers.append(int(groups[1]))
                charges.append(int(groups[2]))
            else:
                raise ValueError(f"String {label} does not match the expected fragment label pattern")
        fragment_types.append(types)
        fragment_numbers.append(numbers)
        fragment_charges.append(charges)

    fragment_type_array = np.array(fragment_types)
    fragment_number_array = np.array(fragment_numbers)
    fragment_charge_array = np.array(fragment_charges)
    mask = np.where((fragment_charge_array > precursor_charges) | (fragment_number_array >= seq_lengths))
    fragment_type_array[mask] = "N"
    fragment_number_array[mask] = 0
    fragment_charge_array[mask] = 0

    return {"type": fragment_type_array, "number": fragment_number_array, "charge": fragment_charge_array}


class SpectralLibrary:
    """
    Main to init a SpectralLibrary obj and go through the steps.

    1- gen_lib
    2- grpc_predict
    3- write output
    """

    library: Spectra
    config: Config
    num_threads: int
    grpc_output: dict

    def __init__(self, search_path: Union[str, Path], out_path: Union[str, Path], config_path: Union[str, Path]):
        """
        Initialize a SpectralLibrary object.

        :param search_path: path to directory containing the msms.txt and raw files
        :param out_path: path to output folder
        :param config_path: path to configuration file
        """
        if isinstance(search_path, str):
            search_path = Path(search_path)
        self.search_path = search_path

        if isinstance(out_path, str):
            out_path = Path(out_path)
        self.out_path = out_path

        self.library = Spectra()

        if isinstance(config_path, str):
            config_path = Path(config_path)
        self.config_path = config_path
        self.config = Config()
        self.config.read(config_path)

        self.out_path.mkdir(exist_ok=True)
        self.results_path = out_path / "results"
        if out_path.exists():
            try:
                self.results_path.mkdir(exist_ok=True)
            except Exception:
                logger.info("In Feature Calculation")
        else:
            logger.info("In Feature Calculation")

    def gen_lib(self):
        """
        Read input csv file and add it to library.

        :raises ValueError: If the value provided for library_input_type in the config file
            is sth. other than "peptides" or "fasta".
        """
        library_input_type = self.config.library_input_type
        if library_input_type == "fasta":
            self.read_fasta()
            library_file = self.out_path / "prosit_input.csv"
        elif library_input_type == "peptides":
            library_file = self.config.library_input
        else:
            raise ValueError(
                f'Library input type {library_input_type} not understood. Can only be "fasta" or "peptides".'
            )
        library_df = csv.read_file(library_file)
        library_df.columns = library_df.columns.str.upper()
        self.library.add_columns(library_df)

    def grpc_predict(self, library: Spectra, alignment: bool = False):
        """
        Use grpc to predict library and add predictions to library.

        :param library: Spectra object with the library
        :param alignment: True if alignment present
        :return: grpc predictions if we are trying to generate spectral library
        """
        triton_client = InferenceServerClient(url=self.config.prediction_server, ssl=self.config.ssl)
        batch_size = 1000

        intensity_outputs = ["intensities", "mz", "annotation"]
        intensity_input_data = {
            "peptide_sequences": (
                library.spectra_data["MODIFIED_SEQUENCE"].to_numpy().reshape(-1, 1).astype(np.object_),
                "BYTES",
            ),
            "collision_energies": (
                library.spectra_data["COLLISION_ENERGY"].to_numpy().reshape(-1, 1).astype(np.float32),
                "FP32",
            ),
            "precursor_charges": (
                library.spectra_data["PRECURSOR_CHARGE"].to_numpy().reshape(-1, 1).astype(np.int32),
                "INT32",
            ),
        }
        intensity_model = self.config.models["intensity"]
        if "tmt" in intensity_model.lower() or "ptm" in intensity_model.lower():
            library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
                lambda x: 2 if x == "HCD" else 1
            )
            intensity_input_data["fragmentation_types"] = (
                library.spectra_data["FRAGMENTATION_GRPC"].to_numpy().reshape(-1, 1).astype(np.float32),
                "FP32",
            )

        intensity_predictions = infer_predictions(
            triton_client,
            model=intensity_model,
            input_data=intensity_input_data,
            outputs=intensity_outputs,
            batch_size=batch_size,
        )
        intensity_predictions["intensities"][np.where(intensity_predictions["intensities"] < 1e-7)] = 0.0

        irt_model = self.config.models["irt"]
        irt_input_data = {"peptide_sequences": intensity_input_data["peptide_sequences"]}
        irt_outputs = ["irt"]
        irt_predictions = infer_predictions(
            triton_client,
            model=irt_model,
            input_data=irt_input_data,
            outputs=irt_outputs,
            batch_size=batch_size,
        )

        if self.config.job_type == "SpectralLibraryGeneration":
            intensity_prediction_dict = {
                "intensity": intensity_predictions["intensities"],
                "fragmentmz": intensity_predictions["mz"],
                "annotation": parse_fragment_labels(
                    intensity_predictions["annotation"],
                    library.spectra_data["PRECURSOR_CHARGE"].to_numpy()[:, None],
                    library.spectra_data["PEPTIDE_LENGTH"].to_numpy()[:, None],
                ),
            }
            output_dict = {intensity_model: intensity_prediction_dict, irt_model: irt_predictions["irt"]}
            return output_dict

        intensities_pred = pd.DataFrame()
        intensities_pred["intensity"] = intensity_predictions["intensities"].tolist()
        library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)

        if alignment:
            return

        library.add_column(irt_predictions["irt"], name="PREDICTED_IRT")

    def read_fasta(self):
        """Read fasta file."""
        cmd = [
            "--fasta",
            f"{self.config.library_input}",
            "--prosit_input",
            f"{self.out_path / 'prosit_input.csv'}",
            "--fragmentation",
            f"{self.config.fragmentation}",
            "--digestion",
            f"{self.config.digestion}",
            "--cleavages",
            f"{self.config.cleavages}",
            "--db",
            f"{self.config.db}",
            "--enzyme",
            f"{self.config.enzyme}",
            "--special-aas",
            f"{self.config.special_aas}",
            "--min-length",
            f"{self.config.min_length}",
            "--max-length",
            f"{self.config.max_length}",
        ]
        digest.main(cmd)
