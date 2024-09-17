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
