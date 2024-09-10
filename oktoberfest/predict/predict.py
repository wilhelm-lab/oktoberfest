import logging
from typing import Dict, List, Optional, Tuple, Union

import anndata
import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics

from ..data.spectra import FragmentType, Spectra
from ..utils import group_iterator
from .koina import Koina

logger = logging.getLogger(__name__)


def predict_intensities(data: anndata.AnnData, chunk_idx: Optional[List[pd.Index]] = None, **kwargs):
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
        >>> pr.predict_intensities(data=library,
        >>>                         model_name="Prosit_2020_intensity_HCD",
        >>>                         server_url="koina.wilhelmlab.org:443",
        >>>                         ssl=True,
        >>>                         targets=["intensities","annotation"])
        >>> print(library.layers["pred_int"])
    """
    if chunk_idx is None:
        intensities = predict_at_once(data=data.obs, **kwargs)
        data.add_intensities(intensities["intensities"], intensities["annotation"], fragment_type=FragmentType.PRED)
    else:
        chunked_intensities = predict_in_chunks(data=data.obs, chunk_idx=chunk_idx, **kwargs)
        data.add_list_of_predicted_intensities(
            chunked_intensities["intensities"], chunked_intensities["annotation"], chunk_idx
        )


def predict_rt(data: anndata.AnnData, **kwargs):
    """
    Retrieve retention time predictions from koina and add them to the provided data object.

    This function takes a dataframe containing information about PSMS and predicts retention time
    using a koina server. The configuration of koina is set using the kwargs.

    :param data: Anndata object containing the data required for prediction and to store the
        predictions in after retrieval from the server.
    :param kwargs: Additional keyword arguments forwarded to Koina::predict

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
        >>> predict_rt(data=library, model_name="Prosit_2019_irt", server_url="koina.wilhelmlab.org:443", ssl=True)
        >>> print(library.obs["PREDICTED_IRT"])
    """
    pred_irts = predict_at_once(data=data.obs, **kwargs)
    data.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")


def predict(
    data: pd.DataFrame, chunk_idx: Optional[List[pd.Index]] = None, **kwargs
) -> Union[Dict[str, List[np.ndarray]], Dict[str, np.ndarray]]:
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

    :return: a dictionary with targets (keys) and predictions (values). If chunk indicies are
        provided, values for each target are a list of numpy array with a length equal to the number
        of chunks provided, else single numpy arrays.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import predict as pr
        >>> import pandas as pd
        >>> # Requiered columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
        >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                         "COLLISION_ENERGY": [30,35],
        >>>                         "PRECURSOR_CHARGE": [1,2],
        >>>                         "FRAGMENTATION": ["HCD","HCD"]})
        >>> predictions = pr.predict(data=meta_df,
        >>>                         model_name="Prosit_2020_intensity_HCD",
        >>>                         server_url="koina.wilhelmlab.org:443",
        >>>                         ssl=True,
        >>>                         targets=["intensities", "annotation"])
        >>> print(predictions)
    """
    if chunk_idx is None:
        return predict_at_once(data, **kwargs)
    return predict_in_chunks(data, chunk_idx, **kwargs)


def predict_at_once(data: pd.DataFrame, **kwargs) -> Dict[str, np.ndarray]:
    """
    Retrieve and return predictions from koina in one go.

    This function takes a dataframe containing information about PSMS and predicts peptide
    properties using a koina server. The configuration of koina is set using the kwargs.
    See the koina predict function for details. TODO, link this properly.

    :param data: Dataframe containing the data for the prediction.
    :param kwargs: Additional keyword arguments forwarded to Koina::predict

    :return: a dictionary with targets (keys) and predictions (values)

    :Example:

    .. code-block:: python

        >>> from oktoberfest import predict as pr
        >>> import pandas as pd
        >>> # Requiered columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE and FRAGMENTATION
        >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                         "COLLISION_ENERGY": [30,35],
        >>>                         "PRECURSOR_CHARGE": [1,2],
        >>>                         "FRAGMENTATION": ["HCD","HCD"]})
        >>> predictions = pr.predict_at_once(data=meta_df,
        >>>                                 model_name="Prosit_2020_intensity_HCD",
        >>>                                 server_url="koina.wilhelmlab.org:443",
        >>>                                 ssl=True,
        >>>                                 targets=["intensities", "annotation"])
        >>> print(predictions)
    """
    predictor = Koina(**kwargs)
    return predictor.predict(data)


def predict_in_chunks(data: pd.DataFrame, chunk_idx: List[pd.Index], **kwargs) -> Dict[str, List[np.ndarray]]:
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

    :Example:

    .. code-block:: python

        >>> from oktoberfest import predict as pr
        >>> from oktoberfest.utils import group_iterator
        >>> # Requiered columns: MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE, FRAGMENTATION and PEPTIDE_LENGTH
        >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                         "COLLISION_ENERGY": [30,35],
        >>>                         "PRECURSOR_CHARGE": [1,2],
        >>>                         "FRAGMENTATION": ["HCD","HCD"],
        >>>                         "PEPTIDE_LENGTH": [8,9]})
        >>> idx = list(group_iterator(df=meta_df, group_by_column="PEPTIDE_LENGTH"))
        >>> predictions = pr.predict_in_chunks(data=meta_df,
        >>>                                 chunk_idx=idx,
        >>>                                 model_name="Prosit_2020_intensity_HCD",
        >>>                                 server_url="koina.wilhelmlab.org:443",
        >>>                                 ssl=True,
        >>>                                 targets=["intensities","annotation"])
        >>> print(predictions)
    """
    predictor = Koina(**kwargs)

    results = []
    for idx in chunk_idx:
        results.append(predictor.predict(data.loc[idx]))
    ret_val = {key: [item[key] for item in results] for key in results[0].keys()}
    return ret_val


def _prepare_alignment_df(library: Spectra, ce_range: Tuple[int, int], group_by_charge: bool = False) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra, and repeating the
    DataFrame for each collision energy (CE) in the given range.

    :param library: the library to be propagated
    :param ce_range: the min and max CE to be propagated for alignment in the dataframe
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :return: a library that is modified according to the description above
    """
    if len(library) < 1000:
        top_n = len(library)
    else:
        top_n = 1000

    if group_by_charge:
        groups = ["RAW_FILE", "PRECURSOR_CHARGE"]
    else:
        groups = ["RAW_FILE"]

    hcd_targets = library.obs.query("(FRAGMENTATION == 'HCD') & ~REVERSE")
    hcd_targets = hcd_targets.sort_values(by="SCORE", ascending=False).groupby(groups)
    top_hcd_targets = hcd_targets.head(top_n)

    alignment_library = library[top_hcd_targets.index]
    alignment_library = Spectra(
        anndata.concat([alignment_library for _ in range(*ce_range)], index_unique="_", keys=range(*ce_range))
    )
    alignment_library.var = library.var
    alignment_library.obs.reset_index(inplace=True)

    alignment_library.obs["ORIG_COLLISION_ENERGY"] = alignment_library.obs["COLLISION_ENERGY"]
    alignment_library.obs["COLLISION_ENERGY"] = np.repeat(range(*ce_range), len(top_hcd_targets))

    return alignment_library


def ce_calibration(
    library: Spectra, ce_range: Tuple[int, int], group_by_charge: bool, model_name: str, **server_kwargs
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
        >>> alignment_library = pr.ce_calibration(library=library,
        >>>                                     ce_range=(15,30),
        >>>                                     group_by_charge=False,
        >>>                                     model_name="Prosit_2020_intensity_HCD",
        >>>                                     server_url="koina.wilhelmlab.org:443",
        >>>                                     ssl=True)
        >>> print(alignment_library)
    """
    alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)

    if "alphapept" in model_name.lower():
        chunk_idx = list(group_iterator(df=alignment_library.obs, group_by_column="PEPTIDE_LENGTH"))
    else:
        chunk_idx = None
    predict_intensities(data=alignment_library, chunk_idx=chunk_idx, model_name=model_name, **server_kwargs)
    _alignment(alignment_library)
    return alignment_library


def _alignment(alignment_library: Spectra):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    """
    pred_intensity = alignment_library.get_matrix(FragmentType.PRED)[0]
    raw_intensity = alignment_library.get_matrix(FragmentType.RAW)[0]
    sm = SimilarityMetrics(pred_intensity, raw_intensity)
    alignment_library.add_column(sm.spectral_angle(raw_intensity, pred_intensity, 0), "SPECTRAL_ANGLE")
