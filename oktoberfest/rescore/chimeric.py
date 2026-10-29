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


def _pick_by_spectral_angle(
    spectrum: np.ndarray,
    pending: np.ndarray,
    obs: np.ndarray,
    mzs: np.ndarray,
    pred: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    peak: np.ndarray,
) -> np.ndarray:
    """
    Pick the candidates by spectral angle on their remaining peaks and remove the peaks of the picked ones (in place).

    :return: whether each candidate lost a peak to a better candidate
    """
    n = len(spectrum)
    pending = pending.copy()
    removed = np.zeros(peak.max() + 1, dtype=bool)
    lost = np.zeros(n, dtype=bool)
    picks = 0
    while (todo := np.flatnonzero(pending)).size:
        picks += 1
        sa = SimilarityMetrics.spectral_angle(obs[todo], pred[todo])
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
    return lost


def _exclusive_rank(
    spectrum: np.ndarray,
    pending: np.ndarray,
    obs: np.ndarray,
    pred: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    peak: np.ndarray,
) -> np.ndarray:
    """
    Rank the candidates of every spectrum by their exclusive spectral angle (0 = best).

    The exclusive spectral angle is the spectral angle on the ions whose observed peak no other candidate of the
    spectrum matched: ions with a shared peak are ignored, predicted but unobserved ions count (NaN counts as 0).
    Ties keep the order of the candidates; candidates of single-candidate spectra get rank 0.
    """
    todo = np.flatnonzero(pending)
    position = np.cumsum(pending) - 1
    # a peak matched by two ions (also of the same candidate) is shared; single-candidate spectra are not ranked
    shared = (np.bincount(peak)[peak] >= 2) & pending[rows]
    exclusive_obs = obs[todo]
    exclusive_pred = pred[todo]
    exclusive_obs[position[rows[shared]], cols[shared]] = 0.0
    exclusive_pred[position[rows[shared]], cols[shared]] = 0.0
    sa = np.nan_to_num(SimilarityMetrics.spectral_angle(exclusive_obs, exclusive_pred))
    order = np.lexsort((-sa, spectrum[todo]))
    first = np.ones(len(order), dtype=bool)
    first[1:] = np.diff(spectrum[todo][order]) != 0
    start = np.maximum.accumulate(np.where(first, np.arange(len(order)), 0))
    rank = np.zeros(len(spectrum), dtype=np.int64)
    rank[todo[order]] = np.arange(len(order)) - start
    return rank


def _pick_by_exclusive_spectral_angle(
    spectrum: np.ndarray,
    pending: np.ndarray,
    obs: np.ndarray,
    mzs: np.ndarray,
    pred: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    peak: np.ndarray,
) -> np.ndarray:
    """
    Pick the candidates by exclusive spectral angle (fixed order) and remove the peaks of the picked ones (in place).

    :return: whether each candidate lost a peak to a better candidate
    """
    rank = _exclusive_rank(spectrum, pending, obs, pred, rows, cols, peak)
    owner = np.full(peak.max() + 1, np.iinfo(np.int64).max)
    np.minimum.at(owner, peak, rank[rows])  # best rank that matched the peak
    gone = (owner[peak] < rank[rows]) & (obs[rows, cols] > EPSILON)
    obs[rows[gone], cols[gone]] = EPSILON
    mzs[rows[gone], cols[gone]] = 0.0
    lost = np.zeros(len(spectrum), dtype=bool)
    lost[rows[gone]] = True
    logger.info(
        f"Subtracted shared peaks of {lost.sum()} candidates in chimeric spectra, picked by exclusive spectral angle"
    )
    return lost


def subtract_shared_peaks(
    metadata: pd.DataFrame,
    raw: sp.csr_matrix,
    mz: sp.csr_matrix,
    pred: sp.csr_matrix,
    order: str = "spectral_angle",
) -> tuple[sp.csr_matrix, sp.csr_matrix, pd.DataFrame]:
    """
    Remove the peaks explained by better candidates of the same spectrum (chimeric spectra).

    Within every spectrum (RAW_FILE, SCAN_NUMBER) with several candidates, the candidate with the highest spectral angle
    is picked, the observed peaks it matched are removed from the spectrum, the spectral angle of the remaining
    candidates is recomputed and the next best one is picked, until all candidates are picked. With order
    "exclusive_spectral_angle", the candidates are picked in a fixed order instead: by their exclusive spectral angle,
    the spectral angle on the ions whose observed peak no other candidate of the spectrum matched (ties in the order of
    the rows). Every candidate keeps the spectrum as it was when it was picked: its matched peaks that a better
    candidate had claimed become "not observed" (intensity EPSILON, no m/z). Observed intensities are relative to the
    most intense annotated peak (MOST_INTESE_PEAK, its absolute intensity); if a better candidate claimed that peak,
    the most intense remaining peak takes its place (intensities rescaled to it, MOST_INTESE_PEAK its absolute
    intensity, NaN if no peak is left), otherwise both stay as they are. All intensity features calculated afterwards
    (spectral angle, observed / predicted counts, fractions and intensities, ...) thus describe the peaks the candidate
    explains on its own. Candidates of single-candidate spectra and the best candidate of every spectrum are unchanged.

    :param metadata: PSM metadata with RAW_FILE, SCAN_NUMBER and optionally MOST_INTESE_PEAK, rows as in the matrices
    :param raw: observed intensities (EPSILON: not observed, 0: invalid ion)
    :param mz: observed m/z of the matched peaks (0 where not observed)
    :param pred: predicted intensities
    :param order: how the candidates are picked: "spectral_angle" (on the remaining peaks) or
        "exclusive_spectral_angle"
    :raises ValueError: if the order is unknown
    :return: observed intensities, m/z and metadata (copy) after the subtraction
    """
    pickers = {"spectral_angle": _pick_by_spectral_angle, "exclusive_spectral_angle": _pick_by_exclusive_spectral_angle}
    if order not in pickers:
        raise ValueError(f"Unknown candidate order {order!r}, must be one of {', '.join(pickers)}.")
    spectrum = metadata.groupby(["RAW_FILE", "SCAN_NUMBER"], sort=False, observed=True).ngroup().to_numpy()
    pending = np.bincount(spectrum)[spectrum] >= 2
    metadata = metadata.copy()
    if not pending.any():
        return raw, mz, metadata

    obs = raw.toarray()
    mzs = mz.toarray()
    pred_dense = pred.toarray()
    rows, cols = np.nonzero(mzs)
    peak = _peak_ids(spectrum, rows, mzs[rows, cols])
    lost = pickers[order](spectrum, pending, obs, mzs, pred_dense, rows, cols, peak)

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
