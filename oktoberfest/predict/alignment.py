import logging

import anndata
import numpy as np
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


def _prepare_alignment_df(
    library: Spectra, ce_range: tuple[int, int], group_by_charge: bool = False, xl: bool = False
) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra for linear and top 20s for cross-linked peptides
    and repeating the DataFrame for each collision energy (CE) in the given range.

    :param library: the library to be propagated
    :param ce_range: the min and max CE to be propagated for alignment in the dataframe
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :param xl: if true, select the top 50 spectra for cross-linked peptide
    :return: a library that is modified according to the description above
    """
    top_n = 1000 if not xl else 20

    if group_by_charge:
        groups = ["RAW_FILE", "PRECURSOR_CHARGE"]
    else:
        groups = ["RAW_FILE"]

    hcd_targets = library.obs.query("(FRAGMENTATION == 'HCD') & ~REVERSE)")
    hcd_targets = hcd_targets[hcd_targets["MODIFIED_SEQUENCE"].str.match(r"^(?:[^U]*(?:\[UNIMOD:737\])?[^U]*)$")]

    hcd_targets = hcd_targets.sort_values(by="SCORE", ascending=False).groupby(groups)
    if len(hcd_targets)<2000:
        top_n= len(hcd_targets)*0.5
    top_hcd_targets = hcd_targets.head(top_n)

    alignment_library = library[top_hcd_targets.index]
    alignment_library = Spectra(
        anndata.concat([alignment_library for _ in range(*ce_range)], index_unique="_", keys=range(*ce_range))
    )
    alignment_library.var = library.var
    alignment_library.obs.reset_index(inplace=True)

    alignment_library.obs["ORIG_COLLISION_ENERGY"] = alignment_library.obs["COLLISION_ENERGY"]
    alignment_library.obs["COLLISION_ENERGY"] = np.repeat(range(*ce_range), len(top_hcd_targets))

    # alignment_library.uns["ion_types"] = library.uns["ion_types"]

    return alignment_library


def _alignment(alignment_library: Spectra, xl: bool = False):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    :param xl: crosslinked or linear peptide
    """
    if xl:
        pred_intensity_a = alignment_library.get_matrix(FragmentType.PRED_A)
        pred_intensity_b = alignment_library.get_matrix(FragmentType.PRED_B)
        raw_intensity_a = alignment_library.get_matrix(FragmentType.RAW_A)
        raw_intensity_b = alignment_library.get_matrix(FragmentType.RAW_B)
        sm_a = SimilarityMetrics(pred_intensity_a, raw_intensity_a)
        sm_b = SimilarityMetrics(pred_intensity_b, raw_intensity_b)
        alignment_library.add_column(sm_a.spectral_angle(raw_intensity_a, pred_intensity_a, 0), "SPECTRAL_ANGLE_A")
        alignment_library.add_column(sm_b.spectral_angle(raw_intensity_b, pred_intensity_b, 0), "SPECTRAL_ANGLE_B")
        alignment_library.add_column(
            (alignment_library.obs["SPECTRAL_ANGLE_A"] + alignment_library.obs["SPECTRAL_ANGLE_B"]) / 2,
            "SPECTRAL_ANGLE",
        )
    else:
        pred_intensity = alignment_library.get_matrix(FragmentType.PRED)
        raw_intensity = alignment_library.get_matrix(FragmentType.RAW)
        sm = SimilarityMetrics(pred_intensity, raw_intensity)
        alignment_library.add_column(sm.spectral_angle(raw_intensity, pred_intensity, 0), "SPECTRAL_ANGLE")
