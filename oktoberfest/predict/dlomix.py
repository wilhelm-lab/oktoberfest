import logging
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from dlomix.constants import PTMS_ALPHABET
from dlomix.data.fragment_ion_intensity import FragmentIonIntensityDataset
from dlomix.losses import masked_pearson_correlation_distance, masked_spectral_distance
from spectrum_fundamentals.constants import FRAGMENTATION_ENCODING
from tensorflow import keras
from tqdm import tqdm

logger = logging.getLogger(__name__)

DUMMY_COLUMN_NAME = "intensities_raw_dummy"


class DLomix:
    """A class for interacting with DLomix models locally for inference."""

    def __init__(self, model_name: str, model_path: Path, output_path: Path):
        self.model_name = model_name
        self.output_path = output_path

        if model_name == "intensity":
            logger.info(f"Loading model weights from {model_path}")
            self.model = keras.models.load_model(
                model_path,
                custom_objects={
                    "masked_spectral_distance": masked_spectral_distance,
                    "masked_pearson_correlation_distance": masked_pearson_correlation_distance,
                },
            )
            self.output_name = "intensities"
        else:
            # TODO implement iRT predictor
            raise NotImplementedError

    def predict(
        self,
        data: Union[dict[str, np.ndarray], pd.DataFrame],
        _async: bool = True,
        debug=False,
    ) -> dict[str, np.ndarray]:

        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame(data)
        processed_data = self.__transform_data(data)

        # TODO check if already exists
        data_path = self.output_path / "dlomix_input.parquet"
        processed_data.to_parquet(data_path)

        # TODO grab & reformat dataset preprocessing logs
        # TODO maximize batch size given memory available
        ds = FragmentIonIntensityDataset(
            test_data_source=str(data_path),
            label_column=DUMMY_COLUMN_NAME,
            model_features=["precursor_charge_onehot", "collision_energy_aligned_normed", "method_nbr"],
            alphabet=PTMS_ALPHABET,
            batch_size=1024,
        )

        # TODO can we extract progress information from the Keras progress bar and nicely integrate it into our progress instead of just hiding it?
        # TODO can we serve this more efficiently? Multi-threading?
        preds = [
            self.model.predict(batch, verbose=0)
            for batch, _ in tqdm(ds.tensor_test_data, desc="Generating predictions")
        ]
        return {self.output_name: np.concatenate(preds)}

    @staticmethod
    def __transform_data(data: pd.DataFrame) -> pd.DataFrame:
        """Transform data into format required by DLomix

        - rename column names from input dataset to corresponding lowercase column names
        - one-hot encode precursor charge
        - scale collision energy
        - integer-encode fragmentation method
        - add dummy label column (required by DLomix dataset classes even for test data)
        """
        return pd.DataFrame(
            {
                "method_nbr": data["FRAGMENTATION"].apply(lambda x: FRAGMENTATION_ENCODING[x]),
                "precursor_charge_onehot": list(np.eye(6)[data["PRECURSOR_CHARGE"].to_numpy() - 1]),
                "collision_energy_aligned_normed": data["COLLISION_ENERGY"] / 100,
                "modified_sequence": data["MODIFIED_SEQUENCE"],
                DUMMY_COLUMN_NAME: list(np.ones((data.shape[0], 30)).astype(np.float64)),
            }
        )
