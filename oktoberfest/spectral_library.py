import logging
import os
from typing import Optional

import numpy as np
import pandas as pd
import tritonclient.grpc as grpcclient
from spectrum_io.file import csv
from spectrum_io.spectral_library import digest

from .constants_dir import CONFIG_PATH
from .data.spectra import FragmentType, Spectra
from .utils.config import Config

logger = logging.getLogger(__name__)


class SpectralLibrary:
    """
    Main to init a SpectralLibrary obj and go through the steps.

    1- gen_lib
    2- grpc_predict
    3- write output
    """

    path: str
    library: Spectra
    config: Config
    config_path: Optional[str]
    num_threads: int
    grpc_output: dict

    def __init__(self, path: str, out_path: str, config_path: Optional[str]):
        """
        Initialize a SpectralLibrary object.

        :param path: path to directory containing the msms.txt and raw files
        :param out_path: path to output folder
        :param config_path: path to configuration file
        """
        self.path = path
        self.library = Spectra()
        self.config_path = config_path
        self.config = Config()
        if config_path:
            self.config.read(config_path)
        else:
            self.config.read(CONFIG_PATH)
        self.results_path = os.path.join(out_path, "results")
        if os.path.isdir(out_path):
            if not os.path.isdir(self.results_path):
                try:
                    os.makedirs(self.results_path)
                except Exception:
                    print("In Feature Calculation")
        else:
            print("In Feature Calculation")

    def gen_lib(self, df_search: Optional[pd.DataFrame] = None):
        """
        Read input csv file and add it to library.

        :param df_search: unused, necessary to ensure same method signature for inheriting function
        """
        if self.config.fasta:
            self.read_fasta()
            library_df = csv.read_file(os.path.join(self.path, "prosit_input.csv"))
        else:
            for file in os.listdir(self.path):
                if file.endswith(".csv"):
                    library_df = csv.read_file(os.path.join(self.path, file))
        library_df.columns = library_df.columns.str.upper()
        self.library.add_columns(library_df)

    def grpc_predict(self, library: Spectra, alignment: bool = False):
        """
        Use grpc to predict library and add predictions to library.

        :param library: Spectra object with the library
        :param alignment: True if alignment present
        :return: grpc predictions if we are trying to generate spectral library
        """
        triton_client = grpcclient.InferenceServerClient(url=self.config.prosit_server)

        result = []
        batch_size = 1000
        num_spec = len(library.spectra_data["MODIFIED_SEQUENCE"])

        # TODO set inputs dynamically based on model
        # CID has no collision energy
        # TMT needs fragmentation
        # ms2pip no ce
        # AlphaPept needs instrument type
        # if tmt_model:
        #     library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
        #         lambda x: 2 if x == "HCD" else 1
        #     )

        for i in range(0, num_spec, batch_size):
            if num_spec < i + batch_size:
                current_batchsize = num_spec - i
            else:
                current_batchsize = batch_size

            inputs = []
            inputs.append(grpcclient.InferInput("peptide_sequences", [current_batchsize, 1], "BYTES"))
            inputs.append(grpcclient.InferInput("collision_energies", [current_batchsize, 1], "FP32"))
            inputs.append(grpcclient.InferInput("precursor_charges", [current_batchsize, 1], "INT32"))
            outputs = [grpcclient.InferRequestedOutput("intensities")]

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

            result.append(
                triton_client.infer(self.config.models["intensity"], inputs=inputs, outputs=outputs).as_numpy(
                    "intensities"
                )
            )

        predictions = np.vstack(result)
        predictions[np.isnan(predictions)] = -1

        intensities_pred = pd.DataFrame()
        intensities_pred["intensity"] = predictions.tolist()
        library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)

        # TODO for spectral library format prediction in dictionary
        # parse new annotation string into three seperate arrays like it was done in prosit_grpc
        # Return only in spectral library generation otherwise add to library
        # if self.config.job_type == "SpectralLibraryGeneration":
        #     return predictions

        if alignment:
            return

        # iRT prediction for Rescoring
        result = []

        for i in range(0, num_spec, batch_size):
            if num_spec < i + batch_size:
                current_batchsize = num_spec - i
            else:
                current_batchsize = batch_size

            inputs = []
            inputs.append(grpcclient.InferInput("peptide_sequences", [current_batchsize, 1], "BYTES"))
            inputs[0].set_data_from_numpy(
                library.spectra_data["MODIFIED_SEQUENCE"]
                .values[i : i + current_batchsize]
                .reshape(-1, 1)
                .astype(np.object_)
            )
            outputs = [grpcclient.InferRequestedOutput("irt")]
            result.append(
                triton_client.infer(self.config.models["irt"], inputs=inputs, outputs=outputs).as_numpy("irt")
            )

        irt_pred = np.vstack(result)
        library.add_column(irt_pred, "PREDICTED_IRT")

    def read_fasta(self):
        """Read fasta file."""
        cmd = [
            "--fasta",
            f"{self.config.fasta}",
            "--prosit_input",
            f"{os.path.join(self.path, 'prosit_input.csv')}",
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
