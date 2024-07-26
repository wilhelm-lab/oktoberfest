from __future__ import annotations

import importlib.util
import logging
from pathlib import Path
from typing import Dict, List

import anndata
import numpy as np
import pandas as pd

from ..data.spectra import FragmentType, Spectra
from ..utils import Config, group_iterator
from .alignment import _alignment, _prepare_alignment_df
from .dlomix import DLomix
from .koina import Koina

logger = logging.getLogger(__name__)


class Predictor:
    """Abstracts common prediction operations away from their actual implementation via the DLomix or Koina interface."""

    def __init__(self, predictor: Koina | DLomix):
        """Initialize from Koina or DLomix instance."""
        self._predictor = predictor

    @classmethod
    def from_config(cls, config: Config, model_type: str, **kwargs) -> Predictor:
        """Load from config object."""
        model_name = config.models[model_type]

        # TODO add documentation for config
        if model_type == "irt" or not config.predict_intensity_locally:
            logger.info(f"Using model {model_name} via Koina")
            koina = Koina(model_name=model_name, server_url=config.prediction_server, ssl=config.ssl, **kwargs)
            return Predictor(koina)

        DLomix.initialize_tensorflow()

        if not importlib.util.find_spec("dlomix"):
            logger.exception(
                ModuleNotFoundError(
                    """Local prediction configured, but the DLomix package could not be found. Please verify that the
                     optional DLomix dependency has been installed."""
                )
            )

        output_folder = config.output / "data/dlomix"

        if model_name == "baseline":
            return Predictor(DLomix(model_type, model_name, output_folder))

        model_path = Path(model_name)
        if model_path.exists():
            logger.info(f"Loading pre-trained PrositIntensityPredictor from {model_path}")
            return Predictor(DLomix(model_type, model_path, output_folder))
        else:
            logger.critical(f"Specified model path {model_name} does not exist, please check")
            logger.exception(FileNotFoundError())

        return

    def predict_intensities(self, data: anndata.AnnData, chunk_idx: List[pd.Index] | None = None, **kwargs):
        """
        Generate intensity predictions and add them to the provided data object.

        This function takes a dataframe containing information about PSMs and predicts intensities. The configuration
        of Koina/DLomix is set using the kwargs. The function either predicts everything at once by concatenating all
        prediction results into single numpy arrays, or returns a list of individual numpy arrays, following the
        indices provided by optionally provided chunks of the dataframe.

        :param data: Anndata object containing the required data for prediction and to store the
            predictions in after retrieval from the server.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional keyword arguments forwarded to Koina/DLomix::predict
        """
        if chunk_idx is None:
            intensities = self.predict_at_once(data=data, **kwargs)
            data.add_intensities(intensities["intensities"], intensities["annotation"], fragment_type=FragmentType.PRED)
        else:
            chunked_intensities = self.predict_in_chunks(data=data.obs, chunk_idx=chunk_idx, **kwargs)
            data.add_list_of_predicted_intensities(
                chunked_intensities["intensities"], chunked_intensities["annotation"], chunk_idx
            )

    def predict_rt(self, data: anndata.AnnData, **kwargs):
        """
        Generate retention time predictions and add them to the provided data object.

        This function takes a dataframe containing information about PSMs and predicts retention times. The
        configuration of Koina/DLomix is set using the kwargs.

        :param data: Anndata object containing the data required for prediction and to store the
            predictions in after retrieval from the server.
        :param kwargs: Additional keyword arguments forwarded to Koina/DLomix::predict
        """
        pred_irts = self.predict_at_once(data=data.obs, **kwargs)
        data.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")

    def predict(
        self, data: pd.DataFrame, chunk_idx: List[pd.Index] | None = None, **kwargs
    ) -> Dict[str, List[np.ndarray]] | Dict[str, np.ndarray]:
        """
        Retrieve and return predictions.

        This function takes a dataframe containing information about PSMs and predicts peptide properties. The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.
        The function either predicts everything at once by concatenating all prediction results
        into single numpy arrays, or returns a list of individual numpy arrays, following the
        indices provided by optionally provided chunks of the dataframe.

        :param data: Dataframe containing the data for the prediction.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and predictions (values). If chunk indices are
            provided, values for each target are a list of numpy array with a length equal to the number
            of chunks provided, else single numpy arrays.
        """
        if chunk_idx is None:
            return self.predict_at_once(data, **kwargs)
        return self.predict_in_chunks(data, chunk_idx, **kwargs)

    def predict_at_once(self, data: Spectra, **kwargs) -> Dict[str, np.ndarray]:
        """
        Retrieve and return predictions in one go.

        This function takes a dataframe containing information about PSMs and predicts peptide properties. The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.

        :param data: Dataframe containing the data for the prediction.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and predictions (values)
        """
        return self._predictor.predict(data, **kwargs)

    def predict_in_chunks(self, data: Spectra, chunk_idx: List[pd.Index], **kwargs) -> Dict[str, List[np.ndarray]]:
        """
        Retrieve and return predictions in chunks.

        This function takes a dataframe containing information about PSMs and predicts peptide properties.The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.

        :param data: Dataframe containing the data for the prediction.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and list of predictions (values) with a length equal
            to the number of chunks provided.
        """
        results = []
        for idx in chunk_idx:
            results.append(self._predictor.predict(data.loc[idx]))
        ret_val = {key: [item[key] for item in results] for key in results[0].keys()}
        return ret_val

    def ce_calibration(
        self, library: Spectra, ce_range: tuple[int, int], group_by_charge: bool, model_name: str, **kwargs
    ) -> Spectra:
        """
        Calculate best collision energy for peptide property predictions.

        The function propagates the provided library object to test NCEs in the given ce range, performs
        intensity prediction for the 1000 highest scoring target PSMs at each NCE and computes the spectral angle
        between predicted and observed intensities before returning the alignment library.

        :param library: spectral library to perform CE calibration on
        :param ce_range: the min and max CE to be tested during calibration
        :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
        :param model_name: The name of the requested prediction model. This is forwarded to the prediction method with
            kwargs and checked here to determine if alphapept is used for further preprocessing.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict
        :return: a spectra object containing the spectral angle for each tested CE
        """
        alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)

        if "alphapept" in model_name.lower():
            chunk_idx = list(group_iterator(df=alignment_library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        self.predict_intensities(data=alignment_library, chunk_idx=chunk_idx, temporary=True, **kwargs)
        _alignment(alignment_library)
        return alignment_library
