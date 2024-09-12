import logging

import anndata
import numpy as np
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


def _prepare_alignment_df(library: Spectra, ce_range: tuple[int, int], group_by_charge: bool = False) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoys and non-HCD-fragmented spectra
    from the input library, selecting the top 1000 (or however many are available if <1000) highest-scoring spectra, and
    repeating the DataFrame for each collision energy (CE) in the given range.

    :param library: the library to be propagated
    :param ce_range: the min and max CE to be propagated for alignment in the dataframe
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :return: a library that is modified according to the description above
    """
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

    alignment_library.uns["ion_types"] = library.uns["ion_types"]

    return alignment_library


def _alignment(alignment_library: Spectra):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    """
    pred_intensity = alignment_library.get_matrix(FragmentType.PRED)
    raw_intensity = alignment_library.get_matrix(FragmentType.RAW)
    sm = SimilarityMetrics(pred_intensity, raw_intensity)
    alignment_library.add_column(sm.spectral_angle(raw_intensity, pred_intensity, 0), "SPECTRAL_ANGLE")
