from __future__ import annotations

import logging
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

from ..data.spectra import FragmentType, Spectra
from ..utils import Config, group_iterator
from .dlomix import DLomix
from .koina import Koina
from .predict import _alignment, _prepare_alignment_df

logger = logging.getLogger(__name__)


class Predictor:
    def __init__(self, predictor: Koina | DLomix):
        self._predictor = predictor

    @classmethod
    def from_config(cls, config: Config, model_name: str, **kwargs) -> Predictor:
        # TODO add documentation for config
        if not config.predict_locally:
            koina = Koina(
                model_name=config.models[model_name], server_url=config.prediction_server, ssl=config.ssl, **kwargs
            )
            return Predictor(koina)
        elif not Path(config.models[model_name]).exists():
            # TODO properly wrap this
            # should this be the default behavior? Or should we throw an exception right away?
            logger.info(f"Path {config.models[model_name]} does not exist, trying to find it on Koina server")
            try:
                return Predictor(
                    Koina(model_name=config.models[model_name], server_url="koina.wilhelmlab.org:443", ssl=True)
                )
                logger.info(f"Successfully using Koina server for model {config.models[model_name]}")
            except Exception as e:
                raise e
        else:
            try:
                import dlomix

                return Predictor(DLomix(model_name, config.models[model_name], config.output / "results", config.inference_batch_size))
            except ImportError as e:
                logger.critical(
                    "DLomix package required for local prediction not found. Please verify that the optional DLomix dependency has been installed."
                )
                raise e

    def predict_intensities(self, data: anndata.AnnData, chunk_idx: list[pd.Index] | None = None, **kwargs):
        """
        Retrieve intensity predictions from koina and add them to the provided data object.

        This function takes a dataframe containing information about PSMS and predicts intensities using
        a koina server. The configuration of koina is set using the kwargs.
        The function either predicts everything at once by concatenating all prediction results
        into single numpy arrays, or returns a list of individual numpy arrays, following the
        indices provided by optionally provided chunks of the dataframe.

        :param data: Anndata object containing the required data for prediction and to store the
            predictions in after retrieval from the server.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional keyword arguments forwarded to Koina::predict
        """
        if chunk_idx is None:
            intensities = self.predict_at_once(data=data.obs, **kwargs)
            data.add_intensities(intensities["intensities"], intensities["annotation"], fragment_type=FragmentType.PRED)
        else:
            chunked_intensities = self.predict_in_chunks(data=data.obs, chunk_idx=chunk_idx, **kwargs)
            data.add_list_of_predicted_intensities(
                chunked_intensities["intensities"], chunked_intensities["annotation"], chunk_idx
            )

    def predict_rt(self, data: anndata.AnnData, **kwargs):
        """
        Retrieve retention time predictions from koina and add them to the provided data object.

        This function takes a dataframe containing information about PSMS and predicts retention time
        using a koina server. The configuration of koina is set using the kwargs.

        :param data: Anndata object containing the data required for prediction and to store the
            predictions in after retrieval from the server.
        :param kwargs: Additional keyword arguments forwarded to Koina::predict
        """
        pred_irts = self.predict_at_once(data=data.obs, **kwargs)
        data.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")

    def predict(
        self, data: pd.DataFrame, chunk_idx: list[pd.Index] | None = None, **kwargs
    ) -> dict[str, list[np.ndarray]] | dict[str, np.ndarray]:
        """
        Retrieve and return predictions from koina.

        This function takes a dataframe containing information about PSMS and predicts peptide
        properties using a koina server. The configuration of koina is set using the kwargs.
        See the koina predict function for details. TODO, link this properly.
        The function either predicts everything at once by concatenating all prediction results
        into single numpy arrays, or returns a list of individual numpy arrays, following the
        indices provided by optionally provided chunks of the dataframe.

        :param data: Dataframe containing the data for the prediction.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional keyword arguments forwarded to Koina::predict

        :return: a dictionary with targets (keys) and predictions (values). If chunk indices are
            provided, values for each target are a list of numpy array with a length equal to the number
            of chunks provided, else single numpy arrays.
        """
        if chunk_idx is None:
            return self.predict_at_once(data, **kwargs)
        return self.predict_in_chunks(data, chunk_idx, **kwargs)

    def predict_at_once(self, data: pd.DataFrame, **kwargs) -> dict[str, np.ndarray]:
        """
        Retrieve and return predictions from koina in one go.

        This function takes a dataframe containing information about PSMS and predicts peptide
        properties using a koina server. The configuration of koina is set using the kwargs.
        See the koina predict function for details. TODO, link this properly.

        :param data: Dataframe containing the data for the prediction.
        :param kwargs: Additonal keyword arguments forwarded to Koina::predict

        :return: a dictionary with targets (keys) and predictions (values)
        """
        return self._predictor.predict(data)

    def predict_in_chunks(self, data: pd.DataFrame, chunk_idx: list[pd.Index], **kwargs) -> dict[str, list[np.ndarray]]:
        """
        Retrieve and return predictions from koina in chunks.

        This function takes a dataframe containing information about PSMS and predicts peptide
        properties using a koina server. The configuration of koina is set using the kwargs.
        See the koina predict function for details. TODO, link this properly.

        :param data: Dataframe containing the data for the prediction.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional keyword arguments forwarded to Koina::predict

        :return: a dictionary with targets (keys) and list of predictions (values) with a length equal
            to the number of chunks provided.
        """
        results = []
        for idx in chunk_idx:
            results.append(self._predictor.predict(data.loc[idx]))
        ret_val = {key: [item[key] for item in results] for key in results[0].keys()}
        return ret_val

    def ce_calibration(
        self, library: Spectra, ce_range: tuple[int, int], group_by_charge: bool, model_name: str, **server_kwargs
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
            server_kwargs and checked here to determine if alphapept is used for further preprocessing.
        :param server_kwargs: Additional parameters that are forwarded to the prediction method
        :return: a spectra object containing the spectral angle for each tested CE
        """
        alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)

        if "alphapept" in model_name.lower():
            chunk_idx = list(group_iterator(df=alignment_library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        self.predict_intensities(data=alignment_library, chunk_idx=chunk_idx, model_name=model_name, **server_kwargs)
        _alignment(alignment_library)
        return alignment_library
