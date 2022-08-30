import logging
import os

import numpy as np
import pandas as pd
from prosit_grpc.predictPROSIT import PROSITpredictor
from prosit_io.file import csv

from .constants import CERTIFICATES, PROSIT_SERVER
from .constants_dir import CONFIG_PATH
from .data.spectra import FragmentType, Spectra
from .utils.config import Config

logger = logging.getLogger(__name__)


class SpectralLibrary:
    """
    main to init a SpectralLibrary obj and go through the steps:
    1- gen_lib
    2- grpc_predict
    3- write output
    """

    path: str
    library: Spectra
    config: dict
    config_path: str
    num_threads: int
    grpc_output: dict

    def __init__(self, path, out_path, config_path=None):
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
                except:
                    print("In Feature Calculation")
        else:
            print("In Feature Calculation")

    def gen_lib(self):
        """
        Read input csv file and add it to library
        """
        if self.config.get_fasta():
            self.read_fasta()
        else:
            print(self.path)
            for file in os.listdir(self.path):
                if file.endswith(".csv"):
                    library_df = csv.read_file(os.path.join(self.path, file))
            library_df.columns = library_df.columns.str.upper()
            self.library.add_columns(library_df)

    def grpc_predict(self, library, alignment=False):
        """
        Use grpc to predict library and add predictions to library
        :return: grpc predictions if we are trying to generate spectral library
        """
        from pathlib import Path

        path = Path(__file__).parent / "certificates/"
        logger.info(path)
        predictor = PROSITpredictor(
            server="proteomicsdb.org:8500",
            path_to_ca_certificate=os.path.join(path, "Proteomicsdb-Prosit-v2.crt"),
            path_to_certificate=os.path.join(path, "oktoberfest-production.crt"),
            path_to_key_certificate=os.path.join(path, "oktoberfest-production.key"),
        )
        models_dict = self.config.get_models()
        models = []
        tmt_model = False
        for key, value in models_dict.items():
            if value:
                if "TMT" in value:
                    tmt_model = True
                models.append(value)
                if alignment:
                    break

        if tmt_model:

            # TODO: find better way instead of hard coded x[12:]
            if self.config.get_tag() == "tmtpro":
                i = 13
            else:
                i = 12
            library.spectra_data["GRPC_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"].apply(lambda x: x[i:])
            library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
                lambda x: 2 if x == "HCD" else 1
            )
            predictions, sequences = predictor.predict(
                sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                collision_energies=library.spectra_data["COLLISION_ENERGY"].values / 100.0,
                fragmentation=library.spectra_data["FRAGMENTATION_GRPC"].values,
                models=models,
                disable_progress_bar=True,
            )
        else:
            library.spectra_data["GRPC_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"]
            try:
                predictions, sequences = predictor.predict(
                    sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                    charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                    collision_energies=library.spectra_data["COLLISION_ENERGY"].values / 100.0,
                    models=models,
                    disable_progress_bar=True,
                )
            except BaseException:
                logger.exception("An exception was thrown!", exc_info=True)
                print(library.spectra_data["GRPC_SEQUENCE"])

        # Return only in spectral library generation otherwise add to library
        if self.config.get_job_type() == "SpectralLibraryGeneration":
            return predictions
        intensities_pred = pd.DataFrame()
        intensities_pred["intensity"] = predictions[models[0]]["intensity"].tolist()

        library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)
        if alignment:
            return
        irt_pred = predictions[models[1]]
        library.add_column(irt_pred, "PREDICTED_IRT")

        if len(models) > 2:
            proteotypicity_pred = predictions[models[2]]
            library.add_column(proteotypicity_pred, "PROTEOTYPICITY")

    def read_fasta(self):
        pass
