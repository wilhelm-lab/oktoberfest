import logging
from enum import Enum
from pathlib import Path
from threading import Thread
from typing import List, Tuple, Type, TypeVar, Union

import numpy as np
import pandas as pd
import scipy
import spectrum_fundamentals.constants as c
from anndata import AnnData
from scipy.sparse import csr_matrix
from spectrum_io.file import hdf5

logger = logging.getLogger(__name__)


SpectraT = TypeVar("SpectraT", bound="Spectra")


class FragmentType(Enum):
    """FragmentType class to enumerate pred, raw, and mz."""

    PRED = 1
    RAW = 2
    MZ = 3


class Spectra:
    """Main to init spectra data."""

    INTENSITY_COLUMN_PREFIX = "INTENSITY_RAW"
    INTENSITY_PRED_PREFIX = "INTENSITY_PRED"
    MZ_COLUMN_PREFIX = "MZ_RAW"
    INTENSITY_PRED_LAYER_NAME = "pred_int"
    INTENSITY_LAYER_NAME = "raw_int"
    MZ_LAYER_NAME = "mz"
    COLUMNS_FRAGMENT_ION = ["Y1+", "Y1++", "Y1+++", "B1+", "B1++", "B1+++"]

    spectra_data: AnnData

    def __init__(self, size):
        """Initialize spectra data as an AnnData object."""
        vars_df = self._gen_vars_df()
        self.spectra_data = AnnData(shape=size, var=vars_df)

    @staticmethod
    def _gen_vars_df() -> pd.DataFrame:
        """
        Creates Annotation dataframe for vars in AnnData object.

        :return: pd.Dataframe of Frgment Annotation
        """
        ion_nums = np.repeat(np.arange(1, 30), 6)
        ion_charge = np.tile([1, 2, 3], 29 * 2)
        temp_cols = []
        for size in range(1, 30):
            for typ in ["Y", "B"]:
                for charge in ["+", "++", "+++"]:
                    temp_cols.append(f"{typ}{size}{charge}")
        ion_types = [frag[0] for frag in temp_cols]
        var_df = pd.DataFrame({"ion": temp_cols, "num": ion_nums, "type": ion_types, "charge": ion_charge})
        var_df = var_df.set_index("ion")
        return var_df

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType) -> List[str]:
        """
        Get column names of the spectra data.

        :param fragment_type: choose predicted, raw, or mz
        :return: A list of column names
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        for i in range(1, 30):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + "_" + column.replace("1", str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type: FragmentType) -> str:
        """
        Resolve prefix given fragment type (1 for pred, 2 for raw, 3 for mz).

        :param fragment_type: choose predicted, raw, or mz
        :return: prefix as string
        """
        logger.debug(fragment_type)
        if fragment_type.value == 1:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type.value == 2:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        elif fragment_type.value == 3:
            prefix = Spectra.MZ_COLUMN_PREFIX
        return prefix

    @staticmethod
    def _resolve_layer_name(fragment_type: FragmentType) -> str:
        """
        Resolve layer name given fragment type (1 for pred, 2 for raw, 3 for mz).

        :param fragment_type: choose predicted, raw, or mz
        :return: layer name as string
        """
        if fragment_type.value == 1:
            layer = Spectra.INTENSITY_PRED_LAYER_NAME
        elif fragment_type.value == 2:
            layer = Spectra.INTENSITY_LAYER_NAME
        elif fragment_type.value == 3:
            layer = Spectra.MZ_LAYER_NAME
        return layer

    def add_column(self, column_data: np.ndarray, name: str) -> None:
        """
        Add column to spectra data.

        :param column_data: data for a column as np
        :param name: name of the column
        """
        column_df = pd.DataFrame({name: list(column_data)})
        self.spectra_data.obs = pd.concat([self.spectra_data.obs.reset_index(drop=True), column_df], axis=1)

    def add_columns(self, columns_data: pd.DataFrame) -> None:
        """
        Assigns columns to the datastructures in AnnData based on what type they are.

        :param columns_data: a pandas data frame to add can be metrics or metadata
        """
        mz_cols = list(filter(lambda c: c.startswith("MZ_RAW"), columns_data.columns))
        raw_int_cols = list(filter(lambda c: c.startswith("INTENSITY_RAW"), columns_data.columns))
        pred_int_cols = list(filter(lambda c: c.startswith("INTENSITY_PRED"), columns_data.columns))
        meta_cols = list(
            filter(lambda c: not (c.startswith("INTENSITY") or c.startswith("MZ_RAW")), columns_data.columns)
        )

        # replaces X and layers (their shape can't be changed)
        if mz_cols:
            self.spectra_data.layers["mz"] = scipy.sparse.csr_matrix(columns_data[mz_cols])
        if pred_int_cols:
            self.spectra_data.layers["pred_int"] = scipy.sparse.csr_matrix(columns_data[pred_int_cols])
        if raw_int_cols:
            self.spectra_data.layers["raw_int"] = scipy.sparse.csr_matrix(columns_data[raw_int_cols])
        if meta_cols:
            self.spectra_data.obs = columns_data[meta_cols]

    def get_meta_data(self) -> pd.DataFrame:
        """Get meta data with intensity, mz and intensity predictions as pd.DataFrame."""
        return self.spectra_data.obs

    def add_matrix_from_hdf5(self, intensity_data: pd.DataFrame, fragment_type: FragmentType) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity sparse matrix
        :param fragment_type: choose predicted, raw, or mz
        """
        # add sparse matrix of intensities to corresponding layer
        layer = self._resolve_layer_name(fragment_type)
        self.spectra_data.layers[layer] = scipy.sparse.csr_matrix(intensity_data)

    def add_matrix(self, intensity_data, fragment_type: FragmentType, annotation=None) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity numpy array to add
        :param fragment_type: choose predicted, raw, or mz
        """
        # Change zeros to epislon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_data[intensity_data == 0] = c.EPSILON
        intensity_data[intensity_data == -1] = 0.0

        intensity_data = csr_matrix(intensity_data)

        layer = self._resolve_layer_name(fragment_type)

        if annotation:
            if self.spectra_data.layers[layer] is None:
                self.spectra_data.layers[layer] = csr_matrix(intensity_data.shape)
            index = [list(self.spectra_data.var_names).index(i) for i in annotation]
            self.spectra_data.layers[layer][:, index] = intensity_data
        else:
            self.spectra_data.layers[layer] = intensity_data

    def get_matrix(self, fragment_type: FragmentType) -> Tuple[csr_matrix, List[str]]:
        """
        Get intensities sparse matrix from AnnData object.

        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """
        prefix = self._resolve_prefix(fragment_type)
        logger.debug(prefix)

        layer = self._resolve_layer_name(fragment_type)
        matrix = self.spectra_data.layers[layer]

        return matrix, self._gen_column_names(fragment_type)

    def write_as_hdf5(self, output_file: Union[str, Path], include_intensities: bool = False) -> Thread:
        """
        Write intensity and mz data as hdf5.

        :param output_file: path to output file
        :return: the thread object from the hdf5 writer for later joining
        """
        logger.info(self.spectra_data.layers)
        logger.info(self.spectra_data.obs.dtypes)
        # self.spectra_data.write(output_file, compression='gzip')

        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY]
        meta_data = self.get_meta_data()
        if include_intensities:
            sparse_matrix_mz, columns_mz = self.get_matrix(FragmentType.MZ)
            sparse_matrix_intensity_raw, columns_intensity = self.get_matrix(FragmentType.RAW)
        else:
            columns_mz = self._gen_column_names(FragmentType.MZ)
            columns_intensity = self._gen_column_names(FragmentType.RAW)
            sparse_matrix_mz = scipy.sparse.csr_matrix(np.zeros((len(meta_data), 174)))
            sparse_matrix_intensity_raw = scipy.sparse.csr_matrix(np.zeros((len(meta_data), 174)))

        data_sets = [meta_data, sparse_matrix_intensity_raw, sparse_matrix_mz]
        column_names = [columns_intensity, columns_mz]

        return hdf5.write_file(data_sets, output_file, data_set_names, column_names)

    def write_pred_as_hdf5(self, output_file: Union[str, Path]) -> Thread:
        """
        Write intensity, mz, and pred data as hdf5.

        :param output_file: path to output file
        :return: the thread object from the hdf5 writer for later joining

        """
        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY, hdf5.INTENSITY_PRED_KEY]

        sparse_matrix_intensity_raw, columns_intensity = self.get_matrix(FragmentType.RAW)
        sparse_matrix_mz, columns_mz = self.get_matrix(FragmentType.MZ)
        sparse_matrix_pred, columns_pred = self.get_matrix(FragmentType.PRED)
        data_sets = [self.get_meta_data(), sparse_matrix_intensity_raw, sparse_matrix_mz, sparse_matrix_pred]
        column_names = [columns_intensity, columns_mz, columns_pred]

        return hdf5.write_file(data_sets, output_file, data_set_names, column_names)

    @classmethod
    def from_hdf5(cls: Type[SpectraT], input_file: Union[str, Path]):
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        """
        input_file = str(input_file)
        data = anndata.read_hdf(input_file,key="raw_mz")
        spectra = cls(data.shape)
        spectra.spectra_data = data
        """
        input_file = str(input_file)
        sparse_raw_intensities = hdf5.read_file(input_file, f"sparse_{hdf5.INTENSITY_RAW_KEY}")
        sparse_raw_mzs = hdf5.read_file(input_file, f"sparse_{hdf5.MZ_RAW_KEY}")
        spectra = cls(sparse_raw_intensities.shape)
        if not sparse_raw_intensities.empty:
            spectra.add_matrix_from_hdf5(sparse_raw_intensities, FragmentType.RAW)
        if not sparse_raw_mzs.empty:
            spectra.add_matrix_from_hdf5(sparse_raw_mzs, FragmentType.MZ)
        spectra.add_columns(hdf5.read_file(input_file, hdf5.META_DATA_KEY))
        return spectra

    @classmethod
    def from_csv(cls: Type[SpectraT], input_file: Union[str, Path]) -> SpectraT:
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        input_file = str(input_file)
        all_columns = pd.read_csv(input_file)
        size = (all_columns.shape[0], 174)
        spectra = cls(size)
        spectra.add_columns(all_columns)

        return spectra

    def convert_to_df(self) -> pd.DataFrame:
        """
        Gives back spectra_data instance as a pandas Dataframe.

        :return: a pandas DataFrame
        """
        df_merged = self.spectra_data.obs
        logger.debug(self.spectra_data.obs.columns)

        if "mz" in list(self.spectra_data.layers):
            mz_cols = pd.DataFrame(self.get_matrix(FragmentType.MZ)[0].toarray())
            mz_cols.columns = self._gen_column_names(FragmentType.MZ)
            df_merged = pd.concat([df_merged,mz_cols], axis=1)
        if "raw_int" in list(self.spectra_data.layers):
            raw_cols = pd.DataFrame(self.get_matrix(FragmentType.RAW)[0].toarray())
            raw_cols.columns = self._gen_column_names(FragmentType.RAW)
            df_merged = pd.concat([df_merged, raw_cols], axis=1)
        if "pred_int" in list(self.spectra_data.layers):
            pred_cols = pd.DataFrame(self.get_matrix(FragmentType.PRED)[0].toarray())
            pred_cols.columns = self._gen_column_names(FragmentType.PRED)
            df_merged = pd.concat([df_merged, pred_cols], axis=1)
        return df_merged
