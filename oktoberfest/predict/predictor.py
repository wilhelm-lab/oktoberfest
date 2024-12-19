from __future__ import annotations

import importlib
import inspect
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, Union

import numpy as np
import pandas as pd

from ..data.spectra import FragmentType, Spectra
from ..utils import Config, group_iterator
from .alignment import _alignment, _prepare_alignment_df
from .koina import Koina
from .utils import ZeroPredictor

logger = logging.getLogger(__name__)

if TYPE_CHECKING or importlib.util.find_spec("dlomix"):
    from .dlomix import DLomix

    PredictionInterface = Union[DLomix, Koina, ZeroPredictor]
elif importlib.util.find_spec("torch"):
    from .torch import TorchModel
    PredictionInterface = Union[TorchModel, Koina, ZeroPredictor]
else:
    PredictionInterface = Union[Koina, ZeroPredictor]
    DLomix = None


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
    ) -> Predictor:
        """Create Koina predictor."""
        return Predictor(
            Koina(
                model_name=model_name,
                server_url=server_url,
                ssl=ssl,
                targets=targets,
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
    def from_torch(
        cls,
        model_path: Union[str, bytes, os.PathLike],
        ion_dict_path: Union[str, bytes, os.PathLike],
        token_dict_path: Union[str, bytes, os.PathLike],
        yaml_dir_path: Union[str, bytes, os.PathLike],
    ) -> Predictor:
        return Predictor(
            TorchModel(
                model_path=model_path,
                ion_dict_path=ion_dict_path,
                token_dict_path=token_dict_path,
                yaml_dir_path=yaml_dir_path,
            ),
            model_name='torch'
        )

    @classmethod
    def from_config(cls, config: Config, model_type: str, **kwargs) -> Predictor:
        """Load from config object."""
        model_name = config.models[model_type]

        if model_type == "irt" and model_name == "zero_irt":
            logger.info("Using zero predictions for iRT")
            return Predictor(ZeroPredictor(), "zero_iRT")

        if model_type == "irt" or not config.predict_intensity_locally:
            logger.info(f"Using model {model_name} via Koina")
            return Predictor.from_koina(
                model_name=model_name, server_url=config.prediction_server, ssl=config.ssl, **kwargs
            )
        
        # TODO actually pass the output folder through kwargs
        output_folder = config.output / "data/dlomix"
        output_folder.mkdir(parents=True, exist_ok=True)
        
        if model_name == "local":
            return Predictor.from_torch(
                model_path=config.models['weights_path'],
                ion_dict_path=config.models['ion_dict_path'],
                token_dict_path=config.models['token_dict_path'],
                yaml_dir_path=config.models['yaml_dir_path'],
            )
        
        if config.download_baseline_intensity_predictor:
            model_path = output_folder / "prosit_baseline_model.keras"
            download = True
        else:
            model_path = Path(model_name)
            download = False
        return Predictor.from_dlomix(
            model_type, model_path, output_folder, config.dlomix_inference_batch_size, download
        )
    
    def _filter_kwargs(self, **kwargs) -> dict[str, Any]:
        """
        Get only arguments accepted by predictor implementation's predict() method from arbitrary set of kwargs.

        :param kwargs: Set of keyword arguments

        :return: Filtered set of keyword arguments
        """
        signature = inspect.signature(self._predictor.predict)
        return {key: value for key, value in kwargs.items() if key in signature.parameters}

    def predict_intensities(self, data: Spectra, chunk_idx: Optional[list[pd.Index]] = None, **kwargs):
        """
        Generate intensity predictions and add them to the provided data object.

        This function takes a Spectra object containing information about PSMs and predicts intensities. The configuration
        of Koina/DLomix is set using the kwargs. The function either predicts everything at once by concatenating all
        prediction results into single numpy arrays, or returns a list of individual numpy arrays, following the
        indices provided by optionally provided chunks of the dataframe.

        :param data: Spectra object containing the required data for prediction and to store the
            predictions in after retrieval from the server.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional keyword arguments forwarded to Koina/DLomix::predict

        :Example:

        .. code-block:: python

            >>> from oktoberfest.data.spectra import Spectra
            >>> from oktoberfest import predict as pr
            >>> # Requiered columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
            >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> library.strings_to_categoricals()
            >>> intensity_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2020_intensity_HCD",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True,
            >>>                         targets=["intensities", "annotation"])
            >>> intensity_predictor.predict_intensities(data=library)
            >>> print(library.layers["pred_int"])
        """
        if chunk_idx is None:
            intensities = self.predict_at_once(data=data, **kwargs)
            data.add_intensities(intensities["intensities"], intensities["annotation"], fragment_type=FragmentType.PRED)
        else:
            chunked_intensities = self.predict_in_chunks(data=data, chunk_idx=chunk_idx, **kwargs)
            data.add_list_of_predicted_intensities(
                chunked_intensities["intensities"], chunked_intensities["annotation"], chunk_idx
            )

    def predict_rt(self, data: Spectra, **kwargs):
        """
        Generate retention time predictions and add them to the provided data object.

        This function takes a Spectra object containing information about PSMs and predicts retention times. The
        configuration of Koina/DLomix is set using the kwargs.

        :param data: Spectra object containing the data required for prediction and to store the
            predictions in after retrieval from the server.
        :param kwargs: Additional keyword arguments forwarded to Koina/DLomix::predict

        :Example:

        .. code-block:: python

            >>> from oktoberfest.data.spectra import Spectra
            >>> from oktoberfest import predict as pr
            >>> import pandas as pd
            >>> # Requiered columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
            >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> library.strings_to_categoricals()
            >>> irt_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2019_irt",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True)
            >>> irt_predictor.predict_rt(data=library)
            >>> print(library.obs["PREDICTED_IRT"])
        """
        pred_irts = self.predict_at_once(data=data, **kwargs)
        data.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")

    def predict_at_once(self, data: Spectra, **kwargs) -> dict[str, np.ndarray]:
        """
        Retrieve and return predictions in one go.

        This function takes a Spectra object containing information about PSMs and predicts peptide properties. The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.

        :param data: Spectra containing the data for the prediction.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and predictions (values)

        :Example:

        .. code-block:: python

            >>> from oktoberfest import predict as pr
            >>> import pandas as pd
            >>> # Required columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
            >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> intensity_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2020_intensity_HCD",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True,
            >>>                         targets=["intensities", "annotation"])
            >>> predictions = intensity_predictor.predict_at_once(data=library)
            >>> print(predictions)
        """
        return self._predictor.predict(data, **self._filter_kwargs(**kwargs))

    def _predict_at_once_df(self, data: pd.DataFrame, **kwargs) -> dict[str, np.ndarray]:
        """
        Retrieve and return predictions in one go.

        This function takes a dataframe containing information about PSMs and predicts peptide properties. The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.

        :param data: dataframe containing the data for the prediction.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and predictions (values)

        :Example:

        .. code-block:: python

            >>> from oktoberfest import predict as pr
            >>> import pandas as pd
            >>> # Required columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
            >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> intensity_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2020_intensity_HCD",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True,
            >>>                         targets=["intensities", "annotation"])
            >>> predictions = intensity_predictor.predict_at_once(data=library)
            >>> print(predictions)
        """
        return self._predictor.predict(data, **self._filter_kwargs(**kwargs))

    def predict_in_chunks(self, data: Spectra, chunk_idx: list[pd.Index], **kwargs) -> dict[str, list[np.ndarray]]:
        """
        Retrieve and return predictions in chunks.

        This function takes a Spectra object containing information about PSMs and predicts peptide properties.The
        configuration of Koina/DLomix is set using the kwargs.
        See the Koina or DLomix predict functions for details. TODO, link this properly.

        :param data: Spectra object containing the data for the prediction.
        :param chunk_idx: The chunked indices of the provided dataframe. This is required in some cases,
            e.g. if padding should be avoided when predicting peptides of different length.
            For alphapept, this is required as padding is only performed within one batch, leading to
            different sizes of arrays between individual prediction batches that cannot be concatenated.
        :param kwargs: Additional parameters that are forwarded to Koina/DLomix::predict

        :return: a dictionary with targets (keys) and list of predictions (values) with a length equal
            to the number of chunks provided.

        :Example:

        .. code-block:: python

            >>> from oktoberfest import predict as pr
            >>> from oktoberfest.utils import group_iterator
            >>> # Required columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE, FRAGMENTATION and PEPTIDE_LENGTH
            >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"],
            >>>                         "PEPTIDE_LENGTH": [8,9]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> idx = list(group_iterator(df=library.obs, group_by_column="PEPTIDE_LENGTH"))
            >>> intensity_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2020_intensity_HCD",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True,
            >>>                         targets=["intensities", "annotation"])
            >>> predictions = intensity_predictor.predict_in_chunks(data=library, chunk_idx=idx)
            >>> print(predictions)
        """
        results = []
        for idx in chunk_idx:
            results.append(self._predictor.predict(data[idx], **self._filter_kwargs(**kwargs)))
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

        :Example:

        .. code-block:: python

            >>> from oktoberfest.data.spectra import FragmentType, Spectra
            >>> from oktoberfest import predict as pr
            >>> import pandas as pd
            >>> import numpy as np
            >>> # Required columns: RAW_FILE, MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE, REVERSE and SCORE
            >>> meta_df = pd.DataFrame({"RAW_FILE": ["File1","File1"],
            >>>                         "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
            >>>                         "COLLISION_ENERGY": [30,35],
            >>>                         "PRECURSOR_CHARGE": [1,2],
            >>>                         "FRAGMENTATION": ["HCD","HCD"],
            >>>                         "REVERSE": [False,False],
            >>>                         "SCORE": [0,0]})
            >>> var = Spectra._gen_vars_df()
            >>> library = Spectra(obs=meta_df, var=var)
            >>> raw_intensities = np.random.rand(2,174)
            >>> annotation = np.array([var.index,var.index])
            >>> library.add_intensities(raw_intensities, annotation, FragmentType.RAW)
            >>> library.strings_to_categoricals()
            >>> intensity_predictor = pr.Predictor.from_koina(
            >>>                         model_name="Prosit_2020_intensity_HCD",
            >>>                         server_url="koina.wilhelmlab.org:443",
            >>>                         ssl=True,
            >>>                         targets=["intensities", "annotation"])
            >>> alignment_library = intensity_predictor.ce_calibration(library=library, ce_range=(15,30), group_by_charge=False)
            >>> print(alignment_library)
        """
        alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)

        if "alphapept" in self.model_name.lower():
            chunk_idx = list(group_iterator(df=alignment_library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        self.predict_intensities(data=alignment_library, chunk_idx=chunk_idx, keep_dataset=False, **kwargs)
        _alignment(alignment_library)
        return alignment_library
