import logging
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Type, TypeVar, Union

import anndata
import numpy as np
import pandas as pd
import scipy
import spectrum_fundamentals.constants as c
from scipy.sparse import csr_matrix

logger = logging.getLogger(__name__)


SpectraT = TypeVar("SpectraT", bound="Spectra")


class FragmentType(Enum):
    """FragmentType class to enumerate pred, raw, and mz."""

    PRED = 1
    PRED_A = 2
    PRED_B = 3
    RAW = 4
    RAW_A = 5
    RAW_B = 6
    MZ = 7
    MZ_A = 8
    MZ_B = 9


class Spectra(anndata.AnnData):
    """Main to init spectra data."""

    INTENSITY_COLUMN_PREFIX = "INTENSITY_RAW"
    INTENSITY_COLUMN_PREFIX_A = "INTENSITY_RAW_A"
    INTENSITY_COLUMN_PREFIX_B = "INTENSITY_RAW_B"
    INTENSITY_PRED_PREFIX = "INTENSITY_PRED"
    INTENSITY_PRED_PREFIX_A = "INTENSITY_PRED_A"
    INTENSITY_PRED_PREFIX_B = "INTENSITY_PRED_B"
    MZ_COLUMN_PREFIX = "MZ_RAW"
    MZ_COLUMN_PREFIX_A = "MZ_RAW_A"
    MZ_COLUMN_PREFIX_B = "MZ_RAW_B"
    INTENSITY_PRED_LAYER_NAME = "pred_int"
    INTENSITY_PRED_LAYER_NAME_A = "pred_int_A"
    INTENSITY_PRED_LAYER_NAME_B = "pred_int_B"
    INTENSITY_LAYER_NAME = "raw_int"
    INTENSITY_LAYER_NAME_A = "raw_int_A"
    INTENSITY_LAYER_NAME_B = "raw_int_B"
    MZ_LAYER_NAME = "mz"
    MZ_LAYER_NAME_A = "mz_A"
    MZ_LAYER_NAME_B = "mz_B"
    COLUMNS_FRAGMENT_ION = ["Y1+", "Y1++", "Y1+++", "B1+", "B1++", "B1+++"]

    @staticmethod
    def _gen_vars_df(xl: bool = False) -> pd.DataFrame:
        """
        Creates Annotation dataframe for vars in AnnData object.

        :return: pd.Dataframe of Frgment Annotation
        """
        if xl:
            max_range = 59
        else:
            max_range = 30
        ion_nums = np.repeat(np.arange(1, max_range), 6)
        ion_charge = np.tile([1, 2, 3], (max_range - 1) * 2)
        temp_cols = []
        for size in range(1, max_range):
            for typ in ["Y", "B"]:
                for charge in ["+", "++", "+++"]:
                    temp_cols.append(f"{typ}{size}{charge}")
        ion_types = [frag[0] for frag in temp_cols]
        var_df = pd.DataFrame({"ion": temp_cols, "num": ion_nums, "type": ion_types, "charge": ion_charge})
        var_df = var_df.set_index("ion")
        return var_df

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType, xl: bool = False) -> List[str]:
        """
        Get column names of the spectra data.

        :param fragment_type: choose predicted, raw, or mz
        :return: A list of column names
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        if xl:
            max_range = 59
        else:
            max_range = 30
        for i in range(1, max_range):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + "_" + column.replace("1", str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type: FragmentType) -> str:
        """
        Resolve prefix given fragment type (1 for pred, 2 for xl_pred_a, 3 for xl_pred_a, 4 for raw, 5 for xl_raw_a, 6 for xl_raw_b, 7 for mz, 8 for xl_mz_a, 9 for xl_mz_b).

        :param fragment_type: choose predicted, raw, or mz
        :return: prefix as string
        """
        logger.debug(fragment_type)
        if fragment_type.value == 1:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type.value == 2:
            prefix = Spectra.INTENSITY_PRED_PREFIX_A
        elif fragment_type.value == 3:
            prefix = Spectra.INTENSITY_PRED_PREFIX_B
        elif fragment_type.value == 4:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        elif fragment_type.value == 5:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX_A
        elif fragment_type.value == 6:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX_B
        elif fragment_type.value == 7:
            prefix = Spectra.MZ_COLUMN_PREFIX
        elif fragment_type.value == 8:
            prefix = Spectra.MZ_COLUMN_PREFIX_A
        else:
            prefix = Spectra.MZ_COLUMN_PREFIX_B
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
            layer = Spectra.INTENSITY_PRED_LAYER_NAME_A
        elif fragment_type.value == 3:
            layer = Spectra.INTENSITY_PRED_LAYER_NAME_B
        elif fragment_type.value == 4:
            layer = Spectra.INTENSITY_COLUMN_PREFIX
        elif fragment_type.value == 5:
            layer = Spectra.INTENSITY_COLUMN_PREFIX_A
        elif fragment_type.value == 6:
            layer = Spectra.INTENSITY_COLUMN_PREFIX_B
        elif fragment_type.value == 7:
            layer = Spectra.MZ_COLUMN_PREFIX
        elif fragment_type.value == 8:
            layer = Spectra.MZ_COLUMN_PREFIX_A
        elif fragment_type.value == 9:
            layer = Spectra.MZ_COLUMN_PREFIX_B
        return layer

    def add_column(self, data: Union[np.ndarray, pd.Series], name: Optional[str] = None) -> None:
        """
        Add column to spectra data.

        :param data: data for a column as np
        :param name: Optional name of the column, required if providing data as a numpy array.
            In case of a pd.Series, providing a name replaces the series' name.

        :raises AssertionError: if data is not 1-dimensional or a column name is missing when
            providing a numpy array.
        :raises TypeError: if the data type is not understood
        """
        if isinstance(data, np.ndarray):
            if name is None:
                raise AssertionError("Missing column name.")
            if data.ndim != 1:
                raise AssertionError("Column data must be supplied as a 1D numpy array.")
            self.obs[name] = data
        elif isinstance(data, pd.Series):
            self.obs[name or data.name] = data
        else:
            raise TypeError(f"Unsupported data type provided: {type(data)}")

    def get_meta_data(self) -> pd.DataFrame:
        """Get meta data with intensity, mz and intensity predictions as pd.DataFrame."""
        return self.obs

    def add_matrix_from_hdf5(self, intensity_data: pd.DataFrame, fragment_type: FragmentType) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity sparse matrix
        :param fragment_type: choose predicted, raw, or mz
        """
        # add sparse matrix of intensities to corresponding layer
        layer = self._resolve_layer_name(fragment_type)
        self.layers[layer] = scipy.sparse.csr_matrix(intensity_data)

    def add_matrix(
        self, intensity_data: np.ndarray, fragment_type: FragmentType, annotation: Optional[np.ndarray] = None
    ) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity numpy array to add with shape (n x m)
        :param fragment_type: choose predicted, raw, or mz
        :param annotation: Optional fragment ion annotations in ProForma notation with shape (n x m)
        """
        # Change zeros to epislon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_data[intensity_data == 0] = c.EPSILON
        intensity_data[intensity_data == -1] = 0.0

        intensity_data = csr_matrix(intensity_data)
        layer = self._resolve_layer_name(fragment_type)

        if annotation:
            if self.layers[layer] is None:
                self.layers[layer] = csr_matrix(intensity_data.shape)
            index = [list(self.var_names).index(i) for i in annotation]
            self.layers[layer][:, index] = intensity_data
        else:
            self.layers[layer] = intensity_data

    def get_matrix(self, fragment_type: FragmentType) -> Tuple[csr_matrix, List[str]]:
        """
        Get intensities sparse matrix from AnnData object.

        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """
        prefix = self._resolve_prefix(fragment_type)
        logger.debug(prefix)

        layer = self._resolve_layer_name(fragment_type)
        matrix = self.layers[layer]

        return matrix, self._gen_column_names(fragment_type)

    def write_as_hdf5(self, output_file: Union[str, Path]):
        """
        Write spectra_data to hdf5 file.

        :param output_file: path to output file
        """
        self.write(output_file, compression="gzip")

    @classmethod
    def from_hdf5(cls: Type[SpectraT], input_file: Union[str, Path]):
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        
        return cls(anndata.read_h5ad(str(input_file)))

    def convert_to_df(self) -> pd.DataFrame:
        """
        Gives back spectra_data instance as a pandas Dataframe.

        :return: a pandas DataFrame
        """
        df_merged = self.obs
        logger.debug(self.obs.columns)

        if "mz" in list(self.layers):
            mz_cols = pd.DataFrame(self.get_matrix(FragmentType.MZ)[0].toarray())
            mz_cols.columns = self._gen_column_names(FragmentType.MZ)
            df_merged = pd.concat([df_merged, mz_cols], axis=1)
        if "raw_int" in list(self.layers):
            raw_cols = pd.DataFrame(self.get_matrix(FragmentType.RAW)[0].toarray())
            raw_cols.columns = self._gen_column_names(FragmentType.RAW)
            df_merged = pd.concat([df_merged, raw_cols], axis=1)
        if "pred_int" in list(self.layers):
            pred_cols = pd.DataFrame(self.get_matrix(FragmentType.PRED)[0].toarray())
            pred_cols.columns = self._gen_column_names(FragmentType.PRED)
            df_merged = pd.concat([df_merged, pred_cols], axis=1)
        return df_merged
