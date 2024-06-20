import logging
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from tensorflow import keras

logger = logging.getLogger(__name__)

class DLomix:
    """A class for interacting with DLomix models locally for inference."""
    def __init__(
        self,
        model_name: str,
        model_path: Path
    ):
        if model_name == "intensity":
            # TODO multi-threading stuff
            logger.info(f"Loading model weights from {model_path}")
            model = keras.models.load_model(model_path)
            """
            # TODO which of these need to be parametrized?
            seq_length = 30
            alphabet = "PTMS_ALPHABET"
            input_mapping = None
            meta_data_keys = None

            self.model = PrositIntensityPredictor(
                seq_length=seq_length,
                alphabet=alphabet,
                use_prosit_ptm_features=False,
                with_termini=False,
                input_keys=input_mapping,
                meta_data_keys=meta_data_keys
            )
            """
            self.model = model
        else:
            # TODO implement iRT predictor
            pass

    def predict(
        self,
        data: Union[dict[str, np.ndarray], pd.DataFrame],
        _async: bool = True,
        debug=False,
    ) -> dict[str, np.ndarray]:
        pass


