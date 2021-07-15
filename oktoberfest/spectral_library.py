import numpy as np
import pandas as pd
import json
import os

from .data.spectra import Spectra
from .data.spectra import FragmentType
from prosit_io.file import csv
from prosit_grpc.predictPROSIT import PROSITpredictor
from .constants import CERTIFICATES, PROSIT_SERVER
from constants_dir import CONFIG_PATH

def read_config():
    print(CONFIG_PATH)
    with open(CONFIG_PATH) as f:
        data = json.load(f)
    return data


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

    def __init__(self, path):
        self.path = path
        self.library = Spectra()
        self.config = read_config()

    def gen_lib(self):
        """
        Read input csv file and add it to library
        """
        if self.config['fileUploads']['fasta']:
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
        predictor = PROSITpredictor(server=PROSIT_SERVER,
                                    path_to_ca_certificate=CERTIFICATES['CA'],
                                    path_to_certificate=CERTIFICATES['USER'],
                                    path_to_key_certificate=CERTIFICATES['KEY'],
                                    keepalive_timeout_ms=10000)

        models_dict = self.config['models']
        models = []
        for key, value in models_dict.items():
            if value:
                models.append(value)

        print(library.spectra_data["PRECURSOR_CHARGE"].values.tolist()[0]>1)
        predictions = predictor.predict(sequences=library.spectra_data["MODIFIED_SEQUENCE"].values,
                                        charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                                        collision_energies=library.spectra_data["COLLISION_ENERGY"].values,
                                        models=models,
                                        disable_progress_bar=True)

    #     return predictions
    #
    # def somefunc():
        #Return only in spectral library generation otherwise add to library
        if self.config['jobType'] == "SpectralLibraryGeneration":
            return predictions

        intensities_pred = pd.DataFrame()
        intensities_pred['intensity'] = predictions[models[0]]["intensity"].tolist()
        library.add_matrix(intensities_pred['intensity'], FragmentType.PRED)
        irt_pred = predictions[models[1]]
        library.add_column(irt_pred, 'IRT')

        if len(models) > 2:
            proteotypicity_pred = predictions[models[2]]
            library.add_column(proteotypicity_pred, 'PROTEOTYPICITY')



    def read_fasta(self):
        pass
