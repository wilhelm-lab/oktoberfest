import numpy as np
import pandas as pd

class SpectralLibrary:
    """
        main to init a SpectralLibrary obj and go through the steps:
        1- gen_lib
        2- grpc_predict
        3- write output
    """
    path: str
    library: pd.DataFrame

    def __init__(self, path):
        self.path = path
        self.library = pd.DataFrame()

    def gen_lib(self):
        """
        Read input csv file and add it to library
        """
        pass

    def grpc_predict(self):
        """
        Use grpc to predict library and add predictions to library
        :return:
        """
        pass