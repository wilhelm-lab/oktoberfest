import logging
import re
from typing import Dict, Tuple

import anndata
import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics

from ..data.spectra import FragmentType, Spectra
from .koina import Koina

logger = logging.getLogger(__name__)


def predict(data: pd.DataFrame, xl : bool = False, **kwargs) -> Dict[str, np.ndarray]:
    """
    Retrieve predictions from koina.

    This function takes a dataframe containing information about PSMS and predicts peptide
    properties using a koina server. The configuration of koina is set using the kwargs.
    See the koina predict function for details. TODO, link this properly.

    :param data: Dataframe containing the data for the prediction.
    :param kwargs: Additional keyword arguments forwarded to Koina::predict

    :return: a dictionary with targets (keys) and predictions (values)
    """
    predictor = Koina(**kwargs)
    if xl:
        results_a, results_b = predictor.predict_xl(data)
        return results_a, results_b 
    else:
        results = predictor.predict(data)
        return results 


def parse_fragment_labels(
    spectra_labels: np.ndarray, precursor_charges: np.ndarray, seq_lengths: np.ndarray
) -> Dict[str, np.ndarray]:
    """Uses regex to parse labels."""
    pattern = rb"([y|b])([0-9]{1,2})\+([1-3])"
    fragment_types = []
    fragment_numbers = []
    fragment_charges = []
    for spectrum_labels in spectra_labels:
        types = []
        numbers = []
        charges = []
        for label in spectrum_labels:
            match = re.match(pattern, label)
            if match:
                groups = match.groups()
                types.append(groups[0].decode())
                numbers.append(int(groups[1]))
                charges.append(int(groups[2]))
            else:
                raise ValueError(f"String {label} does not match the expected fragment label pattern")
        fragment_types.append(types)
        fragment_numbers.append(numbers)
        fragment_charges.append(charges)

    fragment_type_array = np.array(fragment_types)
    fragment_number_array = np.array(fragment_numbers)
    fragment_charge_array = np.array(fragment_charges)
    mask = np.where((fragment_charge_array > precursor_charges) | (fragment_number_array >= seq_lengths))
    fragment_type_array[mask] = "N"
    fragment_number_array[mask] = 0
    fragment_charge_array[mask] = 0

    return {"type": fragment_type_array, "number": fragment_number_array, "charge": fragment_charge_array}


def _prepare_alignment_df(library: Spectra, ce_range: Tuple[int, int], group_by_charge: bool = False, xl : bool = False) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra for linear and top 50s for cross-linked peptides
    and repeating the DataFrame for each collision energy (CE) in the given range.

    :param library: the library to be propagated
    :param ce_range: the min and max CE to be propagated for alignment in the dataframe
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :param xl: if true, select the top 50 spectra for cross-linked peptide
    :return: a library that is modified according to the description above
    """
    top_n = 1000 if not xl else 50
    hcd_targets = library.obs.query("(FRAGMENTATION == 'HCD') & ~REVERSE")
    hcd_targets = hcd_targets.sort_values(by="SCORE", ascending=False)
    if group_by_charge:
        hcd_targets = hcd_targets.groupby("PRECURSOR_CHARGE")
    top_hcd_targets = hcd_targets.head(top_n)

    alignment_library = library[top_hcd_targets.index]
    alignment_library = Spectra(anndata.concat([alignment_library for _ in range(*ce_range)]))
    alignment_library.obs.reset_index(inplace=True)

    alignment_library.obs["ORIG_COLLISION_ENERGY"] = alignment_library.obs["COLLISION_ENERGY"]
    alignment_library.obs["COLLISION_ENERGY"] = np.repeat(range(*ce_range), top_n)

    return alignment_library


def ce_calibration(library: Spectra, ce_range: Tuple[int, int], group_by_charge: bool, xl : bool = False, **server_kwargs ) -> Spectra:
    """
    Calculate best collision energy for peptide property predictions.

    The function propagates the provided library object to test NCEs in the given ce range, performs
    intensity prediction for the 1000 highest scoring target PSMs at each NCE and computes the spectral angle
    between predicted and observed intensities before returning the alignment library.

    :param library: spectral library to perform CE calibration on
    :param ce_range: the min and max CE to be tested during calibration
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :param server_kwargs: Additional parameters that are forwarded to the prediction method
    :param xl: if true, do ce_calibration on cross-linked peptide
    :return: a spectra object containing the spectral angle for each tested CE
    """
    alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge, xl=xl)
    if xl:
        intensities_a, intensities_b = predict(alignment_library.obs, xl = True, **server_kwargs)
        alignment_library.add_matrix(intensities_a["intensities"], FragmentType.PRED_A)
        alignment_library.add_matrix(intensities_b["intensities"], FragmentType.PRED_B)
    else:
        intensities = predict(alignment_library.obs, **server_kwargs)
        alignment_library.add_matrix(intensities["intensities"], FragmentType.PRED)
    
    _alignment(alignment_library, xl=xl)
    return alignment_library


def _alignment(alignment_library: Spectra, xl: bool = False):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    """
    if xl:
        pred_intensity_a = alignment_library.get_matrix(FragmentType.PRED_A)[0]
        pred_intensity_b = alignment_library.get_matrix(FragmentType.PRED_B)[0]
        raw_intensity_a = alignment_library.get_matrix(FragmentType.RAW_A)[0]
        raw_intensity_b = alignment_library.get_matrix(FragmentType.RAW_B)[0]
        sm_a = SimilarityMetrics(pred_intensity_a, raw_intensity_a)
        sm_b = SimilarityMetrics(pred_intensity_b, raw_intensity_b)
        alignment_library.add_column(sm_a.spectral_angle(raw_intensity_a, pred_intensity_a, 0), "SPECTRAL_ANGLE_A")
        alignment_library.add_column(sm_b.spectral_angle(raw_intensity_b, pred_intensity_b, 0), "SPECTRAL_ANGLE_B")
        alignment_library.add_column((alignment_library.obs["SPECTRAL_ANGLE_A"] + alignment_library.obs["SPECTRAL_ANGLE_B"]) / 2, "SPECTRAL_ANGLE")
    else:
        pred_intensity = alignment_library.get_matrix(FragmentType.PRED)[0]
        raw_intensity = alignment_library.get_matrix(FragmentType.RAW)[0]
        sm = SimilarityMetrics(pred_intensity, raw_intensity)
        alignment_library.add_column(sm.spectral_angle(raw_intensity, pred_intensity, 0), "SPECTRAL_ANGLE")


    
    
