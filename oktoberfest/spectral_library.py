import logging
import re
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import tritonclient.grpc as grpcclient
from spectrum_io.file import csv
from spectrum_io.spectral_library import digest

from .constants_dir import CONFIG_PATH
from .data.spectra import FragmentType, Spectra
from .utils.config import Config

logger = logging.getLogger(__name__)


def parse_fragment_labels(spectra_labels: List[np.ndarray]):
    """Uses regex to parse labels."""
    pattern = rb"([y|b])([0-9]{1,2})\+([1-3])"
    annotation = {"type": [], "number": [], "charge": []}
    for spectrum_labels in spectra_labels:
        types = []
        numbers = []
        charges = []
        for label in spectrum_labels:
            match = re.match(pattern, label)
            if match:
                groups = match.groups()
                types.append(groups[0])
                numbers.append(int(groups[1]))
                charges.append(int(groups[2]))
            else:
                raise ValueError(f"String {label} does not match the expected fragment label pattern")
        annotation["type"].append(types)
        annotation["number"].append(numbers)
        annotation["charge"].append(charges)
    for key, value in annotation.items():
        annotation[key] = np.array(value)
    return annotation


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

    def __init__(
        self, search_path: Union[str, Path], out_path: Union[str, Path], config_path: Optional[Union[str, Path]] = None
    ):
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

        if config_path is None:
            config_path = CONFIG_PATH
        if isinstance(config_path, str):
            config_path = Path(config_path)
        self.config_path = config_path
        self.config = Config()
        self.config.read(config_path)

        self.results_path = out_path / "results"
        if out_path.exists():
            try:
                self.results_path.mkdir(exist_ok=True)
            except Exception:
                logger.info("In Feature Calculation")
        else:
            logger.info("In Feature Calculation")

    def gen_lib(self):
        """Read input csv file and add it to library."""
        if self.config.fasta:
            self.read_fasta()
            library_file = self.search_path / "prosit_input.csv"
        else:
            library_file = list(self.search_path.glob("*.csv"))[0]
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
        triton_client = grpcclient.InferenceServerClient(url=self.config.prosit_server, ssl=True)

        batch_size = 1000
        num_spec = len(library.spectra_data["MODIFIED_SEQUENCE"])

        intensity_outputs = ["intensities", "mz", "annotation"]
        irt_outputs = ["irt"]
        predictions_intensity = {output: [] for output in intensity_outputs}
        predictions_irt = {output: [] for output in irt_outputs}

        for i in range(0, num_spec, batch_size):
            if num_spec < i + batch_size:
                current_batchsize = num_spec - i
            else:
                current_batchsize = batch_size

            inputs = []
            inputs.append(grpcclient.InferInput("peptide_sequences", [current_batchsize, 1], "BYTES"))
            inputs.append(grpcclient.InferInput("collision_energies", [current_batchsize, 1], "FP32"))
            inputs.append(grpcclient.InferInput("precursor_charges", [current_batchsize, 1], "INT32"))
            pred_outputs_intensity = []
            for output in intensity_outputs:
                pred_outputs_intensity.append(grpcclient.InferRequestedOutput(output))
            pred_outputs_irt = []
            for output in irt_outputs:
                pred_outputs_irt.append(grpcclient.InferRequestedOutput(output))

            inputs[0].set_data_from_numpy(
                library.spectra_data["MODIFIED_SEQUENCE"]
                .values[i : i + current_batchsize]
                .reshape(-1, 1)
                .astype(np.object_)
            )
            inputs[1].set_data_from_numpy(
                library.spectra_data["COLLISION_ENERGY"]
                .values[i : i + current_batchsize]
                .reshape(-1, 1)
                .astype(np.float32)
            )
            inputs[2].set_data_from_numpy(
                library.spectra_data["PRECURSOR_CHARGE"]
                .values[i : i + current_batchsize]
                .reshape(-1, 1)
                .astype(np.int32)
            )
            prediction_intensity = triton_client.infer(
                self.config.models["intensity"], inputs=inputs, outputs=pred_outputs_intensity
            )
            prediction_irt = triton_client.infer(
                self.config.models["irt"], inputs=[inputs[0]], outputs=pred_outputs_irt
            )

            for output in intensity_outputs:
                predictions_intensity[output].append(prediction_intensity.as_numpy(output))
            for output in irt_outputs:
                predictions_irt[output].append(prediction_irt.as_numpy(output))

        for key, value in predictions_intensity.items():
            predictions_intensity[key] = np.vstack(value)

        for key, value in predictions_irt.items():
            predictions_irt[key] = np.vstack(value)

        if self.config.job_type == "SpectralLibraryGeneration":
            output_dict = {self.config.models["intensity"]: {}, self.config.models["irt"]: predictions_irt["irt"]}
            output_dict[self.config.models["intensity"]]["intensity"] = predictions_intensity["intensities"]
            output_dict[self.config.models["intensity"]]["fragmentmz"] = predictions_intensity["mz"]
            output_dict[self.config.models["intensity"]]["annotation"] = parse_fragment_labels(
                predictions_intensity["annotation"]
            )
            return output_dict

        intensities_pred = pd.DataFrame()
        intensities_pred["intensity"] = predictions_intensity["intensities"].tolist()
        library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)

        if alignment:
            return

        library.add_column(predictions_irt["irt"], name="PREDICTED_IRT")

    def read_fasta(self):
        """Read fasta file."""
        cmd = [
            "--fasta",
            f"{self.config.fasta}",
            "--prosit_input",
            f"{self.search_path / 'prosit_input.csv'}",
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
