from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd
from koinapy.grpc import Koina as _KoinaGRPC

from ..data.spectra import Spectra

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy as np


alternative_column_map = {
    "peptide_sequences": "MODIFIED_SEQUENCE",
    "precursor_charges": "PRECURSOR_CHARGE",
    "collision_energies": "COLLISION_ENERGY",
    "fragmentation_types": "FRAGMENTATION",
    "instrument_types": "INSTRUMENT_TYPES",
}

alternative_column_map_xl = {
    "peptide_sequences_1": "MODIFIED_SEQUENCE_A",
    "peptide_sequences_2": "MODIFIED_SEQUENCE_B",
    "precursor_charges": "PRECURSOR_CHARGE",
    "collision_energies": "COLLISION_ENERGY",
    "fragmentation_types": "FRAGMENTATION",
    "instrument_types": "INSTRUMENT_TYPES",
}

# Create a new mapping with switched keys and values for peptide_sequences_1 and peptide_sequences_2
alternative_column_map_xl_switched = {
    "peptide_sequences_1": "MODIFIED_SEQUENCE_B",
    "peptide_sequences_2": "MODIFIED_SEQUENCE_A",
    **{
        key: value
        for key, value in alternative_column_map_xl.items()
        if key not in ["peptide_sequences_1", "peptide_sequences_2"]
    },
}


class Koina(_KoinaGRPC):
    """Extension of the Koina GRPC class in koinapy, to add required logic for Oktoberfest."""

    def predict(self, data: dict[str, np.ndarray] | pd.DataFrame | Spectra, **kwargs) -> dict[str, np.ndarray]:
        """
        Perform inference on the given data using the Koina model.

        This method allows you to perform inference on the provided input data using the configured Koina model. You can
        choose to perform inference asynchronously (in parallel) or sequentially, depending on the value of the '_async'
        parameter. If asynchronous inference is selected, the method will return when all inference tasks are complete.
        Note: Ensure that the model and server are properly configured and that the input data matches the model's
        input requirements.

        :param data: A dictionary or dataframe containing input data for inference. For the dictionary, keys are input names,
            and values are numpy arrays. In case of a dataframe, the input fields for the requested model must be present
            in the column names.
        :param kwargs: Additional params that are forwarded to super().predict
        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.

        Example::
            model = Koina("Prosit_2019_intensity")
            input_data = {
                "peptide_sequences": np.array(["PEPTIDEK" for _ in range(size)]),
                "precursor_charges": np.array([2 for _ in range(size)]),
                "collision_energies": np.array([20 for _ in range(size)]),
                "fragmentation_types": np.array(["HCD" for _ in range(size)]),
                "instrument_types": np.array(["QE" for _ in range(size)])
            }
            predictions = model.predict(input_data)
        """
        if isinstance(data, Spectra):
            data = data.obs
        if isinstance(data, pd.DataFrame):
            data = {
                input_field: data[[alternative_column_map[input_field]]].to_numpy()
                for input_field in self.model_inputs.keys()
            }
        return super().predict(inputs=data, **kwargs)

    def predict_xl(
        self, data: dict[str, np.ndarray] | pd.DataFrame | Spectra, **kwargs
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """
        Perform inference on the xl data using the Koina model.

        This method allows you to perform inference on the provided input data using the configured Koina model. You can
        choose to perform inference asynchronously (in parallel) or sequentially, depending on the value of the '_async'
        parameter. If asynchronous inference is selected, the method will return when all inference tasks are complete.
        Note: Ensure that the model and server are properly configured and that the input data matches the model's
        nput requirements.

        :param data: A dictionary or dataframe containing input data for inference. For the dictionary, keys are input names,
            and values are numpy arrays. In case of a dataframe, the input fields for the requested model must be present
            in the column names.
        :param kwargs: Additional params that are forwarded to super().predict
        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.
        :raises ValueError: If `data` is not of type `Spectra`, `pd.DataFrame`, or a dictionary.

        Example::
            model = Koina("Prosit_XL_CMS2_intensity")
            input_data = {
                "peptide_sequences_1": np.array(["PEPTIDEK" for _ in range(size)]),
                "peptide_sequences_2": np.array(["PEPTIDEK" for _ in range(size)]),
                "precursor_charges": np.array([2 for _ in range(size)]),
                "collision_energies": np.array([20 for _ in range(size)]),
                "fragmentation_types": np.array(["HCD" for _ in range(size)]),
                "instrument_types": np.array(["QE" for _ in range(size)])
            }
            predictions = model.predict(input_data)
        """
        if isinstance(data, Spectra):
            data = data.obs
        if isinstance(data, pd.DataFrame):
            data = {
                input_field: data[[alternative_column_map_xl[input_field]]].to_numpy()
                for input_field in self.model_inputs.keys()
            }
            prediction_ab = super().predict(inputs=data, debug=True, **kwargs)
            temp_field = data["peptide_sequences_1"].copy()
            data["peptide_sequences_1"] = data["peptide_sequences_2"]
            data["peptide_sequences_2"] = temp_field
            prediction_ba = super().predict(inputs=data, debug=True, **kwargs)

            return prediction_ab, prediction_ba

        raise ValueError("Input data must be of type Spectra, pd.DataFrame, or a dictionary.")
