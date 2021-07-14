import os

import numpy as np
import pandas as pd

from .spectral_library import SpectralLibrary
from .data.spectra import Spectra


class CeCalibration(SpectralLibrary):
    """
        main to init a CeCalibrarion obj and go through the steps:
        1- gen_lib
        2- allign_ce
        3- get_best_ce
        4- write output
    """
    msms_path: str
    raw_path: str
    best_ce: float

    def __init__(self, search_path, raw_path):
        self.search_path = search_path
        self.raw_path = raw_path
        self.library = Spectra()

    def gen_lib(self):
        """
        Read input msms and raw and add it to library
        """
        if os.path.basename(self.search_path) == "msms.txt":
            print(True)

    def allign_ce(self):
        """
        Edit library to try different ranges of ce 15-50.
        then predict with the new library.
        Check https://gitlab.lrz.de/proteomics/prosit_tools/oktoberfest/-/blob/develop/oktoberfest/ce_calibration/grpc_alignment.py
        """
        pass

    def get_best_ce(self):
        """
        Get aligned ce for this lib.
        """
        pass

if __name__ == "main":
    pass
