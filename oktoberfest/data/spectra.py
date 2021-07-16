import pandas as pd
import scipy
from scipy.sparse import coo_matrix, spmatrix
from enum import Enum
import numpy as np


class FragmentType(Enum):
    PRED = 1
    RAW = 2
    MZ = 3


class Spectra:
    NO_OF_FRAGMENTS = 174
    INTENSITY_COLUMN_PREFIX = 'INTENSITY_RAW'
    INTENSITY_PRED_PREFIX = 'INTENSITY_PRED'
    MZ_COLUMN_PREFIX = 'MZ_RAW'
    EPSILON = 1e-7
    COLUMNS_FRAGMENT_ION = ['Y1+', 'Y1++', 'Y1+++', 'B1+', 'B1++', 'B1+++']

    spectra_data: pd.DataFrame

    def __init__(self):
        self.spectra_data = pd.DataFrame()

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType):
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        for i in range(1, 30):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + '_' + column.replace('1', str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type):
        if fragment_type is FragmentType.PRED:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type is FragmentType.RAW:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        else:
            prefix = Spectra.MZ_COLUMN_PREFIX
        return prefix

    def add_column(self, column_data: np, name: str):
        self.spectra_data[name] = column_data

    def add_columns(self, columns_data: pd.DataFrame):
        """
        Add columns to spectra data.
        :param columns_data: a pandas data frame to add can be metrics or metadata.
        """
        # Check if columns already exist
        self.spectra_data = pd.concat([self.spectra_data, columns_data], axis=1)

    def get_meta_data(self):
        columns = self.spectra_data.columns
        meta_data_columns = list(filter(lambda x: not x.startswith(
            tuple([Spectra.INTENSITY_COLUMN_PREFIX, Spectra.MZ_COLUMN_PREFIX, Spectra.INTENSITY_PRED_PREFIX])), columns))

        return self.spectra_data[meta_data_columns]

    def add_matrix(self, intensity_data, fragment_type=FragmentType.PRED):
        """
        concat intensity df as a sparse matrix to our data
        :param intensity_data: Intensity numpy array to add
        :param fragment_type: Choose type of fragments predicted or raw
        """
        intensity_df = intensity_data.explode()

        # reshape based on the number of fragments
        intensity_array = intensity_df.values.astype(np.float16).reshape(-1, Spectra.NO_OF_FRAGMENTS)

        # Change zeros to epislon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_array[intensity_array == 0] = Spectra.EPSILON
        intensity_array[intensity_array == -1] = 0

        # generate column names and build dataframe from sparse matrix
        intensity_df = pd.DataFrame.sparse.from_spmatrix(coo_matrix(intensity_array)).astype(np.float16)
        columns = self._gen_column_names(fragment_type)
        intensity_df.columns = columns
        self.add_columns(intensity_df)

    def get_matrix(self, fragment_type=FragmentType.PRED, return_column_names=False) -> spmatrix:
        """
        Get intensities sparse matrix from dataframe.
        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """

        prefix = Spectra._resolve_prefix(fragment_type)
        columns_to_select = list(filter(lambda c: c.startswith(prefix), self.spectra_data.columns))
        if return_column_names:
            return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values), columns_to_select
        # Check if conversion is low change to coo then csr from coo
        return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values)
