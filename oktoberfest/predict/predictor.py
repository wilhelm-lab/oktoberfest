from __future__ import annotations

import importlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
import pandas as pd

from ..data.spectra import FragmentType, Spectra
from ..utils import Config, group_iterator
from .alignment import _alignment, _prepare_alignment_df
from .koina import Koina

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .dlomix import DLomix

    PredictionInterface = Optional[Union[DLomix, Koina]]

else:
    if importlib.util.find_spec("dlomix"):
        from .dlomix import DLomix

        PredictionInterface = Optional[Union[DLomix, Koina]]
    else:
        PredictionInterface = Koina


class Predictor:
    """Abstracts common prediction operations away from their actual implementation via the DLomix or Koina interface."""

    def __init__(self, predictor: PredictionInterface, model_name: str):
        """Initialize from Koina or DLomix instance."""
        self._predictor = predictor
        self.model_name = model_name

    @classmethod
    def from_koina(
        cls,
        model_name: str,
        server_url: str = "koina.wilhelmlab.org:443",
        ssl: bool = True,
        targets: Optional[list[str]] = None,
        disable_progress_bar: bool = False,
    ) -> Predictor:
        """Create Koina predictor."""
        return Predictor(
            Koina(
                model_name=model_name,
                server_url=server_url,
                ssl=ssl,
                targets=targets,
                disable_progress_bar=disable_progress_bar,
            ),
            model_name=model_name,
        )

    @classmethod
    def from_dlomix(
        cls, model_type: str, model_path: Path, output_path: Path, batch_size: int, download: bool = False
    ) -> Predictor:
        """Create DLomix predictor."""
        return Predictor(
            DLomix(
                model_type=model_type,
                model_path=model_path,
                output_path=output_path,
                batch_size=batch_size,
                download=download,
            ),
            model_name=model_path.stem,
        )

    @classmethod
    def from_config(cls, config: Config, model_type: str, **kwargs) -> Predictor:
        """Load from config object."""
        model_name = config.models[model_type]

        if model_type == "irt" and model_name == "zero_irt":
            logger.info("Using zero predictions for iRT")
            return Predictor(None, "zero_iRT")

        if model_type == "irt" or not config.predict_intensity_locally:
            logger.info(f"Using model {model_name} via Koina")
            return Predictor.from_koina(
                model_name=model_name, server_url=config.prediction_server, ssl=config.ssl, **kwargs
            )

        # TODO actually pass the output folder through kwargs
        output_folder = config.output / "data/dlomix"
        output_folder.mkdir(parents=True, exist_ok=True)

        if config.download_baseline_intensity_predictor:
            model_path = output_folder / "prosit_baseline_model.keras"
            download = True
        else:
            model_path = Path(model_name)
            download = False
        return Predictor.from_dlomix(
            model_type, model_path, output_folder, config.dlomix_inference_batch_size, download
        )

    def predict_intensities(self, data: Spectra, chunk_idx: Optional[list[pd.Index]] = None, **kwargs):
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

    def predict_rt(self, data: Spectra, **kwargs):
        """
        Generate retention time predictions and add them to the provided data object.

        This function takes a dataframe containing information about PSMs and predicts retention times. The
        configuration of Koina/DLomix is set using the kwargs.

        :param data: Anndata object containing the data required for prediction and to store the
            predictions in after retrieval from the server.
        :param kwargs: Additional keyword arguments forwarded to Koina/DLomix::predict
        """
        if self._predictor is None:
            data.add_column(data.obs["RETENTION_TIME"], name="PREDICTED_IRT")
        else:
            pred_irts = self.predict_at_once(data=data.obs, **kwargs)
            data.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")

    def predict_at_once(self, data: Spectra, **kwargs) -> dict[str, np.ndarray]:
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

    def predict_in_chunks(self, data: Spectra, chunk_idx: list[pd.Index], **kwargs) -> dict[str, list[np.ndarray]]:
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
            results.append(self._predictor.predict(data[idx], **kwargs))
        ret_val = {key: [item[key] for item in results] for key in results[0].keys()}
        return ret_val

    def ce_calibration(self, library: Spectra, ce_range: tuple[int, int], group_by_charge: bool, **kwargs) -> Spectra:
        """
        Calculate best collision energy for peptide property predictions.

        The function propagates the provided library object to test NCEs in the given ce range, performs
        intensity prediction for the 1000 highest scoring target PSMs at each NCE and computes the spectral angle
        between predicted and observed intensities before returning the alignment library.

        :param library: spectral library to perform CE calibration on
        :param ce_range: the min and max CE to be tested during calibration
        :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict
        :return: a spectra object containing the spectral angle for each tested CE
        """
        alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)

        if "alphapept" in self.model_name.lower():
            chunk_idx = list(group_iterator(df=alignment_library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        self.predict_intensities(data=alignment_library, chunk_idx=chunk_idx, keep_dataset=False, **kwargs)
        _alignment(alignment_library)
        return alignment_library
