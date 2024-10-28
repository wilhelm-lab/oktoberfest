import numpy as np

from oktoberfest.data import Spectra


class ZeroPredictor:
    """Predictor implementation that returns unmodified peptide properties."""

    def __init__(self):
        """Nothing to do here."""
        pass

    def predict(self, data: Spectra) -> dict[str, np.ndarray]:
        """Return unmodified peptide properties.

        :param data: spectral library whose features will be returned

        :return: a dictionary containing original peptide features in the same format as predictions

        """
        # TODO implement for intensity
        return {"irt": data.obs["RETENTION_TIME"]}
    
    def predict_xl(self, data: Spectra) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """
        Perform inference on the xl data using the using a zero predictor.

        This is currently not implemented.

        :param data: spectral library to predict features for
        :raises NotImplementedError: Always.
        """
        raise NotImplementedError("This method is not implemeted")
