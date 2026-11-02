from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import scipy.sparse as sp
from spectrum_fundamentals.constants import EPSILON
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics

logger = logging.getLogger(__name__)

__all__ = ["subtract_shared_peaks"]


def _peak_ids(spectrum: np.ndarray, rows: np.ndarray, mzs: np.ndarray) -> np.ndarray:
    """Number the observed peaks: the same m/z in the same spectrum is the same peak."""
    order = np.lexsort((mzs, spectrum[rows]))
    new = np.ones(len(order), dtype=bool)
    new[1:] = (np.diff(spectrum[rows][order]) != 0) | (np.diff(mzs[order]) != 0)
    peak = np.empty(len(order), dtype=np.int64)
    peak[order] = np.cumsum(new) - 1
    return peak


def subtract_shared_peaks(
    metadata: pd.DataFrame,
    raw: sp.csr_matrix,
    mz: sp.csr_matrix,
    pred: sp.csr_matrix,
) -> tuple[sp.csr_matrix, sp.csr_matrix, pd.DataFrame]:
    """
    Remove the peaks explained by better candidates of the same spectrum (chimeric spectra).

    Within every spectrum (RAW_FILE, SCAN_NUMBER) with several candidates, the candidate with the highest spectral angle
    is picked, the observed peaks it matched are removed from the spectrum, the spectral angle of the remaining
    candidates is recomputed and the next best one is picked, until all candidates are picked. Every candidate keeps
    the spectrum as it was when it was picked: its matched peaks that a better candidate had claimed become "not
    observed" (intensity EPSILON, no m/z). Observed intensities are relative to the most intense annotated peak
    (MOST_INTESE_PEAK, its absolute intensity); if a better candidate claimed that peak, the most intense remaining
    peak takes its place (intensities rescaled to it, MOST_INTESE_PEAK its absolute intensity, NaN if no peak is left),
    otherwise both stay as they are. All intensity features calculated afterwards (spectral angle, observed /
    predicted counts, fractions and intensities, ...) thus describe the peaks the candidate explains on its own.
    Candidates of single-candidate spectra and the best candidate of every spectrum are unchanged.

    :param metadata: PSM metadata with RAW_FILE, SCAN_NUMBER and optionally MOST_INTESE_PEAK, rows as in the matrices
    :param raw: observed intensities (EPSILON: not observed, 0: invalid ion)
    :param mz: observed m/z of the matched peaks (0 where not observed)
    :param pred: predicted intensities
    :return: observed intensities, m/z and metadata (copy) after the subtraction
    """
    spectrum = metadata.groupby(["RAW_FILE", "SCAN_NUMBER"], sort=False, observed=True).ngroup().to_numpy()
    n = len(spectrum)
    pending = np.bincount(spectrum)[spectrum] >= 2
    metadata = metadata.copy()
    if not pending.any():
        return raw, mz, metadata

    obs = raw.toarray()
    mzs = mz.toarray()
    pred_dense = pred.toarray()
    rows, cols = np.nonzero(mzs)
    peak = _peak_ids(spectrum, rows, mzs[rows, cols])
    removed = np.zeros(peak.max() + 1, dtype=bool)
    lost = np.zeros(n, dtype=bool)

    picks = 0
    while (todo := np.flatnonzero(pending)).size:
        picks += 1
        sa = SimilarityMetrics.spectral_angle(obs[todo], pred_dense[todo])
        order = np.lexsort((-sa, spectrum[todo]))
        first = np.ones(len(order), dtype=bool)
        first[1:] = np.diff(spectrum[todo][order]) != 0
        chosen = todo[order[first]]
        pending[chosen] = False
        is_chosen = np.zeros(n, dtype=bool)
        is_chosen[chosen] = True
        removed[peak[is_chosen[rows]]] = True
        gone = removed[peak] & pending[rows] & (obs[rows, cols] > EPSILON)
        obs[rows[gone], cols[gone]] = EPSILON
        mzs[rows[gone], cols[gone]] = 0.0
        lost[rows[gone]] = True
    logger.info(f"Subtracted shared peaks of {lost.sum()} candidates in chimeric spectra ({picks} picks)")

    # the most intense remaining peak replaces a claimed most intense peak
    observed = obs[lost] > EPSILON
    remaining = np.where(observed, obs[lost], 0).max(axis=1)
    scale = np.where(remaining > 0, remaining, 1.0)
    obs[lost] = np.where(observed, obs[lost] / scale[:, None], obs[lost])
    if "MOST_INTESE_PEAK" in metadata:
        most_intense = metadata["MOST_INTESE_PEAK"].to_numpy(dtype=float)
        most_intense[lost] = np.where(remaining > 0, most_intense[lost] * remaining, np.nan)
        metadata["MOST_INTESE_PEAK"] = most_intense
    return sp.csr_matrix(obs), sp.csr_matrix(mzs), metadata
