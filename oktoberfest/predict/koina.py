import logging
from typing import Dict, Union

import numpy as np
import pandas as pd
from koinapy.grpc import Koina as _KoinaGRPC

from ..data.spectra import Spectra

logger = logging.getLogger(__name__)


alternative_column_map = {
    "peptide_sequences": "MODIFIED_SEQUENCE",
    "precursor_charges": "PRECURSOR_CHARGE",
    "collision_energies": "COLLISION_ENERGY",
    "fragmentation_types": "FRAGMENTATION",
    "instrument_types": "INSTRUMENT_TYPES",
}


class Koina(_KoinaGRPC):
    """Extension of the Koina GRPC class in koinapy, to add required logic for Oktoberfest."""

    def predict(self, data: Union[Dict[str, np.ndarray], pd.DataFrame, Spectra], **kwargs) -> Dict[str, np.ndarray]:
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
                input_field: data[alternative_column_map[input_field]].to_numpy()
                for input_field in self.model_inputs.keys()
            }
        return super().predict(inputs=data, **kwargs)
