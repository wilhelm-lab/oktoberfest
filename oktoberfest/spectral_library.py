import numpy as np
import pandas as pd
import os

from .data.spectra import Spectra
from .data.spectra import FragmentType
from prosit_io.file import csv
from prosit_grpc.predictPROSIT import PROSITpredictor
from .constants import CERTIFICATES, PROSIT_SERVER
from .constants_dir import CONFIG_PATH
from .utils.config import Config


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

    def __init__(self, path, config_path=None):
        self.path = path
        self.library = Spectra()
        self.config_path = config_path
        self.config = Config()
        if config_path:
            self.config.read(config_path)
        else:
            self.config.read(CONFIG_PATH)

    def gen_lib(self):
        """
        Read input csv file and add it to library
        """
        if self.config.get_fasta():
            self.read_fasta()
        else:
            for root, dirs, files in os.walk(self.path):
                for file in files:
                    if file.endswith(".csv"):
                         library_df = csv.read_file(self.path + file)
            library_df.columns = library_df.columns.str.upper()
            self.library.add_columns(library_df)

    def grpc_predict(self, library):
        """
        Use grpc to predict library and add predictions to library
        :return: grpc predictions if we are trying to generate spectral library
        """
        predictor = PROSITpredictor(server=self.config.get_prosit_server())
                                    #path_to_ca_certificate=CERTIFICATES['CA'],
                                    #path_to_certificate=CERTIFICATES['USER'],
                                    #path_to_key_certificate=CERTIFICATES['KEY'],
                                    #keepalive_timeout_ms=10000)

        models_dict = self.config.get_models()
        models = []
        tmt_model = False
        for key, value in models_dict.items():
            if value:
                if 'TMT' in value:
                    tmt_model = True
                models.append(value)

        if tmt_model:
            # TODO: find better way instead of hard coded x[12:]
            library.spectra_data['GRPC_SEQUENCE'] = library.spectra_data['MODIFIED_SEQUENCE'].apply(
                lambda x: x[12:])
            predictions,sequences = predictor.predict(sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                                            charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                                            collision_energies=library.spectra_data["COLLISION_ENERGY"].values/100.0,
                                            fragmentation= library.spectra_data["FRAGMENTATION"].values,
                                            models=models,
                                            disable_progress_bar=True)
        else:
            library.spectra_data['GRPC_SEQUENCE'] = library.spectra_data['MODIFIED_SEQUENCE']
            predictions = predictor.predict(sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                                            charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                                            collision_energies=library.spectra_data["COLLISION_ENERGY"].values/100.0,
                                            models=models,
                                            disable_progress_bar=True)

        #Return only in spectral library generation otherwise add to library
        if self.config.get_job_type() == "SpectralLibraryGeneration":
            return predictions

        intensities_pred = pd.DataFrame()
        intensities_pred['intensity'] = predictions[models[0]]["intensity"].tolist()

        library.add_matrix(intensities_pred['intensity'], FragmentType.PRED)
        irt_pred = predictions[models[1]]
        library.add_column(irt_pred, 'PREDICTED_IRT')

        if len(models) > 2:
            proteotypicity_pred = predictions[models[2]]
            library.add_column(proteotypicity_pred, 'PROTEOTYPICITY')

    def read_fasta(self):
        pass

