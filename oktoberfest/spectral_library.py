import logging
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from prosit_grpc.predictPROSIT import PROSITpredictor
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
        cert_path = Path(__file__).parent / "certificates"
        logger.info(cert_path)

        predictor = PROSITpredictor(
            server=self.config.prosit_server,
            path_to_ca_certificate=cert_path / "Proteomicsdb-Prosit-v2.crt",
            path_to_certificate=cert_path / "oktoberfest-production.crt",
            path_to_key_certificate=cert_path / "oktoberfest-production.key",
        )

        models_dict = self.config.models
        models = []
        tmt_model = False
        for _, value in models_dict.items():
            if not value:
                continue
            tmt_model = True if "TMT" in value else tmt_model
            models.append(value)
            if alignment:
                break

        if tmt_model:
            library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
                lambda x: 2 if x == "HCD" else 1
            )

        library.spectra_data.loc[:, "GRPC_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"]
        try:
            predictions = predictor.predict(
                sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                collision_energies=library.spectra_data["COLLISION_ENERGY"].values / 100.0,
                fragmentation=library.spectra_data["FRAGMENTATION_GRPC"].values if tmt_model else None,
                models=models,
                disable_progress_bar=True,
            )
            # Return only in spectral library generation otherwise add to library
            if self.config.job_type == "SpectralLibraryGeneration":
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
        except Exception as e:
            logger.error(e)

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
