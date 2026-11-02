from __future__ import annotations

import itertools
import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .rescore import rescore_with_percolator

logger = logging.getLogger(__name__)

__all__ = [
    "select_decoy_candidates",
    "localization_probabilities",
    "estimate_localization_flr",
    "hybrid_localization_flr",
    "calibrate_localization_score",
    "localization_groups",
    "neighbor_classes",
    "decoy_proxy_fractions",
    "add_adjacent_decoy_candidates",
    "localize_ptms",
    "sdi_features",
    "ambiguity_flags",
]

# features used to score the localization candidates of a spectrum against each other
LOCALIZATION_FEATURES = [
    "fraction_observed",
    "fraction_observed_mod",
    "sum_observed_and_predicted",
    "observed_intensity",
    "spectral_angle",
    "fraction_observed_and_predicted",
    "observed_intensity_mod",
    "sum_observed_and_predicted_mod",
    "fraction_observed_but_not_predicted",
    "NL_FRACTION",
    "NL_NOT_POSSIBLE",
    "observed_intensity_mod_relative",
    "sum_observed_and_predicted_mod_relative",
]

# features that can only support a localization candidate: their percolator weights must not be negative
LOCALIZATION_POSITIVE_FEATURES = [
    "spectral_angle",
    "fraction_observed_mod",
    "observed_intensity_mod",
    "sum_observed_and_predicted_mod",
    "observed_intensity_mod_relative",
    "sum_observed_and_predicted_mod_relative",
    "NL_FRACTION",
]

# residues decoy candidates can be predicted on (decoy_proxy "search"), per unimod id; all sites for other PTMs.
# Phosphorylated decoys are predicted as pS or pT, which lose phosphoric acid like the annotated decoys
DECOY_PROXY_RESIDUES = {21: ["S", "T"]}

# distances of localization candidates to the called sites are capped at this value
MAX_SITE_DISTANCE = 5

# neutral loss features of decoy candidates without the neutral loss on the decoy site (plain decoys), used to score
# a plain copy of every decoy candidate if the decoys are annotated with the neutral loss (decoy_neutral_loss)
PLAIN_DECOY_FEATURES = {"NL_FRACTION": "NL_FRACTION_PLAIN", "NL_NOT_POSSIBLE": "NL_NOT_POSSIBLE_PLAIN"}
PLAIN_DECOY_SUFFIX = "-plain_decoy"
# residue groups of the called sites whose FLR is estimated from the plain decoys: Y does not lose phosphoric acid,
# so plain decoys represent its wrong localizations; S and T (and mixed) use the decoys with the neutral loss
PLAIN_DECOY_GROUPS = ("Y",)

# site-determining ion (SDI) features of every localization candidate against its best competitor (sdi_features), from
# the stored predictions and annotation of the spectra; the first two can only support a candidate
SDI_FEATURES = ["sdi_pred_obs_frac", "sdi_pred_missing", "sdi_none"]
SDI_POSITIVE_FEATURES = ["sdi_pred_obs_frac", "sdi_pred_missing"]
# annotated intensities (raw_int) of the stored annotation: 0 = impossible ion, possible but unobserved ions are marked
# with a small value (1e-7), observed ions are above this threshold
OBSERVED_INTENSITY = 1e-5

# sharpness of the model based error probabilities of the "hybrid" FLR estimator, in units of the slope of the PEP fit
# of the training step: the high one is used for FLRs above about 1-2%, the low one at the low end. Calls without a
# possible site next to the called site (neighbor_classes "no neighbor") use the sharper pair: their candidates differ
# by many site-determining ions, so the scores separate them better than the training slope suggests
MODEL_SHARPNESS_HIGH = 3.0
MODEL_SHARPNESS_LOW = 2.0
MODEL_SHARPNESS_HIGH_NO_NEIGHBOR = 5.0
MODEL_SHARPNESS_LOW_NO_NEIGHBOR = 2.75
MODEL_ERROR_SHARPNESS = {
    "model_error_high": MODEL_SHARPNESS_HIGH,
    "model_error_low": MODEL_SHARPNESS_LOW,
    "model_error_high_no_neighbor": MODEL_SHARPNESS_HIGH_NO_NEIGHBOR,
    "model_error_low_no_neighbor": MODEL_SHARPNESS_LOW_NO_NEIGHBOR,
}


def _ptm_ids(ptm_sites: dict[int, list[str]]) -> str:
    """Regex alternation of the unimod ids of the localized PTMs, e.g. 21|7."""
    return "|".join(str(unimod_id) for unimod_id in ptm_sites)


def _decoy_site_pattern(ptm_sites: dict[int, list[str]], decoy_residues: Optional[list[str]]) -> Optional[str]:
    """Regex matching PTMs localized on decoy residues, e.g. A[UNIMOD:21]; None if no decoy residues are used."""
    if not decoy_residues:
        return None
    return rf"[{''.join(decoy_residues)}]\[UNIMOD:(?:{_ptm_ids(ptm_sites)})\]"


def decoy_proxy_fractions(search_results: pd.DataFrame, ptm_sites: dict[int, list[str]]) -> dict[str, dict[str, float]]:
    """
    Get the frequencies of the residues decoy candidates are predicted on, from the localizations of the search engine.

    Every target PSM is counted once with the sites localized by the search engine (SEARCH_MODIFIED_SEQUENCE), only
    residues in DECOY_PROXY_RESIDUES (all sites of other PTMs).

    :param search_results: search results with the search engine localization in SEARCH_MODIFIED_SEQUENCE
    :param ptm_sites: the localized PTMs {unimod id: [residues]}
    :return: {unimod id: {residue: fraction}}, empty if the search engine localization is not available
    """
    if "SEARCH_MODIFIED_SEQUENCE" not in search_results.columns:
        return {}
    psms = search_results[~search_results["REVERSE"].astype(bool)] if "REVERSE" in search_results else search_results
    psms = psms.drop_duplicates(["RAW_FILE", "SCAN_NUMBER", "SEARCH_MODIFIED_SEQUENCE"])
    fractions = {}
    for unimod_id, sites in ptm_sites.items():
        residues = [aa for aa in DECOY_PROXY_RESIDUES.get(unimod_id, sites) if aa in sites]
        counts = {
            aa: int(psms["SEARCH_MODIFIED_SEQUENCE"].str.count(rf"{aa}\[UNIMOD:{unimod_id}\]").sum()) for aa in residues
        }
        total = sum(counts.values())
        if total:
            fractions[str(unimod_id)] = {aa: count / total for aa, count in counts.items()}
    logger.info(f"Decoy candidates are predicted on the residues with the frequencies {fractions}")
    return fractions


_TOKEN = re.compile(r"([A-Z])((?:\[UNIMOD:\d+\])*)")


def _decoy_sequences(
    sequence: str,
    ptm_id: int,
    sites: list[str],
    priority_residues: list[str],
    adjacent_residues: list[str],
    max_distance: int,
) -> list[str]:
    """
    Get the decoy localization candidates of a localized sequence (see add_adjacent_decoy_candidates).

    :param sequence: one localization candidate of the PSM, in internal format
    :param ptm_id: unimod id of the localized PTM
    :param sites: residues the PTM can be located on
    :param priority_residues: decoy residues by priority; all positions of the first one in the peptide are decoy sites
    :param adjacent_residues: residues that are decoy sites within max_distance of a possible site
    :param max_distance: largest distance of an adjacent decoy site to a possible site
    :return: decoy candidates: one PTM on a decoy site, the others on all combinations of possible sites
    """
    ptm = f"[UNIMOD:{ptm_id}]"
    prefix = re.match(r"^(?:\[UNIMOD:\d+\])*-?", sequence).group(0)
    tokens = [[aa, mods.replace(ptm, "")] for aa, mods in _TOKEN.findall(sequence[len(prefix) :])]
    n_ptm = sequence.count(ptm)
    if n_ptm == 0:
        return []
    residues = [aa for aa, _ in tokens]
    possible = [i for i, aa in enumerate(residues) if aa in sites]
    decoy_positions = set()
    priority = next((aa for aa in priority_residues if aa in residues), None)
    if priority is not None:
        decoy_positions.update(i for i, aa in enumerate(residues) if aa == priority)
    # decoy sites next to the possible sites (not on the first residue)
    for i in possible:
        for j in range(max(1, i - max_distance), min(len(residues), i + max_distance + 1)):
            if residues[j] in adjacent_residues:
                decoy_positions.add(j)
    decoys = []
    for decoy in sorted(decoy_positions):
        for real in itertools.combinations(possible, n_ptm - 1):
            modified = [list(token) for token in tokens]
            for i in (decoy, *real):
                modified[i][1] += ptm
            decoys.append(prefix + "".join(aa + mods for aa, mods in modified))
    return decoys


def add_adjacent_decoy_candidates(
    search_results: pd.DataFrame,
    ptm_sites: dict[int, list[str]],
    priority_residues: list[str],
    adjacent_residues: list[str],
    max_distance: int = 2,
) -> pd.DataFrame:
    """
    Add decoy localization candidates next to the possible sites.

    Wrong localizations are mostly on neighboring sites, which decoys on a fixed residue rarely represent. For every PSM, decoy sites are all positions of the first decoy residue by priority in the peptide (e.g. A, as in
    select_decoy_candidates) and every residue in adjacent_residues within max_distance of a possible site (not the
    first residue). Every decoy candidate carries one PTM on a decoy site and the other PTMs on a combination of
    possible sites. Only phosphorylation is supported; other PTMs get no decoys.

    :param search_results: search results in internal format with the localization candidates on the possible sites
    :param ptm_sites: the localized PTMs {unimod id: [residues]}
    :param priority_residues: decoy residues by priority
    :param adjacent_residues: residues used as decoy sites next to possible sites
    :param max_distance: largest distance of an adjacent decoy site to a possible site
    :return: search results with the decoy candidates added
    """
    if 21 not in ptm_sites:
        return search_results
    psm = ["RAW_FILE", "SCAN_NUMBER", "SEQUENCE", "PRECURSOR_CHARGE"]
    templates = search_results[search_results["MODIFIED_SEQUENCE"].str.contains(r"\[UNIMOD:21\]", regex=True)]
    templates = templates.drop_duplicates(psm)
    unique = templates["MODIFIED_SEQUENCE"].unique()
    decoys = {
        sequence: _decoy_sequences(sequence, 21, ptm_sites[21], priority_residues, adjacent_residues, max_distance)
        for sequence in unique
    }
    added = templates.assign(MODIFIED_SEQUENCE=templates["MODIFIED_SEQUENCE"].map(decoys))
    added = added.explode("MODIFIED_SEQUENCE").dropna(subset=["MODIFIED_SEQUENCE"])
    logger.info(
        f"Added {len(added)} decoy localization candidates ({len(added) / max(len(templates), 1):.1f} per PSM) on the "
        f"residues {priority_residues} (by priority) and {adjacent_residues} within {max_distance} of possible sites"
    )
    return pd.concat([search_results, added], ignore_index=True)


def select_decoy_candidates(
    search_results: pd.DataFrame, ptm_sites: dict[int, list[str]], decoy_residues: list[str]
) -> pd.DataFrame:
    """
    Keep the decoy localization candidates of one decoy residue per peptide with at most one decoy site.

    The decoy residues are used in the given order of priority: every peptide uses the first decoy residue it contains
    (e.g. A, or G if the peptide has no A). Candidates with more than one decoy site or with a decoy site on another
    residue are removed. Candidates without a decoy site are kept.

    :param search_results: search results in internal format with all localization candidates
    :param ptm_sites: the localized PTMs {unimod id: [residues]}
    :param decoy_residues: decoy residues in order of priority
    :return: search results with the selected decoy candidates
    """
    sequences = search_results["SEQUENCE"].astype(str)
    chosen = pd.Series(np.nan, index=search_results.index, dtype=object)
    for residue in reversed(decoy_residues):
        chosen[sequences.str.contains(residue, regex=False)] = residue
    decoy_sites = search_results["MODIFIED_SEQUENCE"].str.findall(
        rf"([{''.join(decoy_residues)}])\[UNIMOD:(?:{_ptm_ids(ptm_sites)})\]"
    )
    keep = [
        len(sites) == 0 or (len(sites) == 1 and sites[0] == residue)
        for sites, residue in zip(decoy_sites, chosen, strict=True)
    ]
    selected = search_results[keep]
    logger.info(
        f"Using decoy residues {decoy_residues} (by priority): kept {int((decoy_sites[keep].map(len) > 0).sum())} "
        f"decoy localization candidates, removed {len(search_results) - len(selected)}"
    )
    return selected


def _average_fold_weights(weights_file: Path, output_file: Path):
    """
    Average the weights of all cross validation bins of a percolator weights file into a single set of weights.

    Percolator --init-weights only uses one set of weights. Averaging the bins makes the scoring independent of the
    random assignment of PSMs to the cross validation bins.

    :param weights_file: percolator weights file (per bin: feature names, normalized weights, raw weights)
    :param output_file: path of the weights file with the averaged weights
    """
    lines = weights_file.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")]
    names = data[0]
    normalized = np.mean([np.array(data[i].split("\t"), dtype=float) for i in range(1, len(data), 3)], axis=0)
    raw = np.mean([np.array(data[i].split("\t"), dtype=float) for i in range(2, len(data), 3)], axis=0)
    output_file.write_text(
        "\n".join(comments + [names, "\t".join(map(str, normalized)), "\t".join(map(str, raw))]) + "\n"
    )


def _constrain_weights(weights_file: Path, output_file: Path, positive_features: Optional[list[str]] = None):
    """
    Set negative percolator weights of LOCALIZATION_POSITIVE_FEATURES (or positive_features) to 0.

    The weights file contains, for every cross validation bin, the feature names, the normalized and the raw weights.

    :param weights_file: percolator weights file
    :param output_file: path of the constrained weights file
    :param positive_features: features whose weights must not be negative, LOCALIZATION_POSITIVE_FEATURES if None
    """
    if positive_features is None:
        positive_features = LOCALIZATION_POSITIVE_FEATURES
    lines = weights_file.read_text().splitlines()
    names = None
    constrained = []
    for line in lines:
        fields = line.split("\t")
        if line.startswith("#"):
            constrained.append(line)
        elif names is None or fields == names:
            names = fields
            constrained.append(line)
        else:
            constrained.append(
                "\t".join(
                    str(max(0.0, float(value))) if name in positive_features else value
                    for name, value in zip(names, fields, strict=True)
                )
            )
    output_file.write_text("\n".join(constrained) + "\n")


def _scan_numbers(psms: pd.DataFrame) -> list[int]:
    """Get scan numbers from the PSMId (<raw file>-<scan number>-...) of percolator results (raw file as str)."""
    return [
        int(psm_id[len(str(raw_file)) :].split("-")[1])
        for psm_id, raw_file in zip(psms["PSMId"], psms["filename"], strict=True)
    ]


def _site_positions(peptides: pd.Series, ptm_pattern: str) -> list[frozenset]:
    """Get the positions (1-based) of the residues carrying a PTM matched by ptm_pattern."""
    site = re.compile(ptm_pattern)
    # remove flanking residues (K.PEPTIDE.R) and N-terminal modifications
    core = peptides.str.replace(r"^[^.]*\.|\.[^.]*$", "", regex=True).str.replace(
        r"^(?:\[UNIMOD:\d+\])+-", "", regex=True
    )
    return [
        frozenset(
            i + 1
            for i, (_, mods) in enumerate(re.findall(r"([A-Z])((?:\[UNIMOD:\d+\])*)", peptide))
            if site.search(mods)
        )
        for peptide in core
    ]


def _site_distance(sites: frozenset, called: frozenset) -> int:
    """
    Get the distance of a localization candidate to the called sites.

    The distance is the largest distance of a moved site to its closest replaced site, capped at MAX_SITE_DISTANCE.
    """
    if not isinstance(called, frozenset):
        return 0
    moved, replaced = sites - called, called - sites
    if not moved or not replaced:
        return 0
    return min(MAX_SITE_DISTANCE, max(min(abs(p - q) for q in replaced) for p in moved))


def _decoy_weights(psms: pd.DataFrame, decoys: pd.Series, ptm_pattern: str, weighting: str = "count") -> pd.Series:
    """
    Expected number of false localizations per spectrum in which a decoy candidate wins.

    "count" (default): every wrong real candidate is as likely to beat the called candidate as a decoy candidate, so
    a decoy win stands for (number of real candidates - 1) / (number of decoy candidates) false localizations.
    "distance": wrong localizations are most often next to the true site, decoy sites can be anywhere in the peptide.
    The chance q(d) that a candidate at distance d to the called sites scores above the called candidate is estimated
    from the decoy candidates of all spectra, and a decoy win stands for sum(q(d) of the wrong real candidates) /
    sum(q(d) of the decoy candidates) false localizations.

    :param psms: all scored candidates (real and decoy) with filename, ScanNr, peptide and score_prob
    :param decoys: mask of the decoy candidates
    :param ptm_pattern: regex matching the localized PTMs
    :param weighting: "count" or "distance"
    :return: per spectrum (filename, ScanNr) the expected false localizations of a decoy win
    """
    spectrum = ["filename", "ScanNr"]
    if weighting == "count":
        keys = [psms[column] for column in spectrum]
        n_decoy = decoys.groupby(keys).sum()
        n_real = (~decoys).groupby(keys).sum()
        n_decoy = n_decoy[n_decoy > 0]
        return ((n_real.reindex(n_decoy.index) - 1).clip(lower=0) / n_decoy).astype(float).rename("decoy_weight")
    psms = psms.assign(sites=_site_positions(psms["peptide"], ptm_pattern))
    # the called candidate of every spectrum: the best real candidate
    real = psms[~decoys]
    called = real.loc[real.groupby(spectrum)["score_prob"].idxmax(), spectrum + ["sites", "score_prob"]]
    psms = psms.merge(
        called.rename(columns={"sites": "called_sites", "score_prob": "called_score"}), on=spectrum, how="left"
    ).set_index(psms.index)
    distance = pd.Series(
        [_site_distance(s, c) for s, c in zip(psms["sites"], psms["called_sites"], strict=True)], index=psms.index
    )
    beats_call = psms["score_prob"] > psms["called_score"]
    # q(d): how often a decoy candidate at distance d to the called sites scores above the called candidate
    q_decoys = beats_call[decoys].groupby(distance[decoys]).mean()
    q = distance.map(q_decoys).fillna(beats_call[decoys].mean()).clip(lower=1e-6)
    wrong_real = ~decoys & (distance > 0)
    wrong = q[wrong_real].groupby([psms.loc[wrong_real, column] for column in spectrum]).sum()
    decoy = q[decoys].groupby([psms.loc[decoys, column] for column in spectrum]).sum()
    return (wrong.reindex(decoy.index).fillna(0) / decoy).rename("decoy_weight")


def _score_slope(psms_file: Path) -> float:
    """
    Get the slope b of the log odds of a percolator result: log((1 - PEP) / PEP) = a + b * score.

    The training step (best candidate per spectrum) has a well behaved PEP fit, unlike the scoring of all candidates
    of all spectra, where percolator's PEP estimate can collapse to a step function for large inputs.

    :param psms_file: percolator psms file with score and posterior_error_prob
    :return: the slope b, 1.0 if the PEPs cannot be fitted
    """
    psms = pd.read_csv(psms_file, sep="\t", usecols=["score", "posterior_error_prob"])
    pep = psms["posterior_error_prob"].to_numpy()
    valid = (pep > 1e-6) & (pep < 1 - 1e-6)
    if valid.sum() < 10 or np.ptp(psms["score"].to_numpy()[valid]) == 0:
        logger.warning(f"Cannot fit the PEPs of {psms_file}, using the scores as log odds for the localization.")
        return 1.0
    return float(np.polyfit(psms["score"].to_numpy()[valid], np.log((1 - pep[valid]) / pep[valid]), 1)[0])


def localization_probabilities(
    candidates: pd.DataFrame,
    train_dir: Path,
    score_dir: Path,
    decoy_pattern: Optional[str] = None,
    ptm_pattern: str = r"\[UNIMOD:21\]",
    probability: str = "pep",
    decoy_weighting: str = "count",
    extra_features: Optional[list[str]] = None,
    extra_positive_features: Optional[list[str]] = None,
    keep_scores: bool = False,
    num_threads: int = 3,
) -> pd.DataFrame:
    """
    Score all localization candidates of each spectrum and compute localization probabilities.

    First, percolator is trained on localization specific features of all candidates, where the candidates of a
    spectrum compete against each other (train_dir). Second, all candidates are scored with the trained weights
    (percolator --static), with a unique ExpMass per candidate so that percolator keeps all candidates of a spectrum
    instead of only the best one (score_dir). The localization probability of a candidate is its probability
    normalized over all candidates of the same spectrum. With probability "pep" (default), it is the inverse posterior
    error probability of the scoring step: loc_prob = (1 / PEP) / sum(1 / PEP). With probability "score", the PEPs of
    the scoring step are not used, since percolator's PEP estimate on all candidates of all spectra can collapse to a
    step function for large inputs; instead the log odds a + b * score of the PEP fit of the training step are used:
    loc_prob = exp(b * score) / sum(exp(b * score)).

    Decoy candidates (PTM on a decoy residue, matched by decoy_pattern) are passive: they are not used for training
    and not for the localization probabilities, so the localization is the same as without decoys. They are only
    scored to flag spectra in which a decoy candidate would have won (decoy_win), used to estimate the FLR. If the
    candidates have plain decoy features (PLAIN_DECOY_FEATURES, decoys annotated with the neutral loss), every decoy
    candidate is also scored as a plain decoy, without the neutral loss on the decoy site (decoy_win_plain).

    :param candidates: rescore features of all localization candidates (one row per candidate)
    :param train_dir: directory for the percolator training on the localization features
    :param score_dir: directory for scoring all candidates with the trained weights
    :param decoy_pattern: regex matching decoy candidates in the Peptide column, None if no decoys are used
    :param ptm_pattern: regex matching the localized PTMs, used to match decoys and wrong candidates by site distance
    :param extra_features: features used in addition to LOCALIZATION_FEATURES (e.g. SDI_FEATURES)
    :param extra_positive_features: extra features whose weights must not be negative
    :param keep_scores: also return the peptide and the percolator score of every candidate
    :param num_threads: number of threads used by percolator
    :return: one row per (non decoy) candidate with PSMId, filename, ScanNr, loc_prob, loc_prob_delta (difference
        between the best and second best candidate of the spectrum), possible_sites (number of candidates of the
        spectrum), decoy_win (a decoy candidate scores higher than all candidates of the spectrum) and decoy_weight,
        and decoy_win_plain and decoy_weight_plain for the plain decoys (same as decoy_win / decoy_weight without
        them), and peptide and score with keep_scores
    """
    features = LOCALIZATION_FEATURES + list(extra_features or [])
    columns = (
        ["SpecId", "Label", "ScanNr", "filename", "ExpMass"]
        + [feature for feature in features if feature in candidates.columns]
        + ["Peptide", "Proteins"]
    )
    df_candidates = candidates[columns].copy()

    if decoy_pattern:
        n_decoy_sites = df_candidates["Peptide"].str.count(decoy_pattern)
        # at most one decoy site per candidate
        df_candidates = df_candidates[n_decoy_sites <= 1]
        is_decoy = n_decoy_sites[n_decoy_sites <= 1] == 1
    else:
        is_decoy = pd.Series(False, index=df_candidates.index)
    use_plain = decoy_pattern is not None and all(
        plain in candidates.columns for plain in PLAIN_DECOY_FEATURES.values()
    )
    if use_plain:
        # plain copies of the decoy candidates: the neutral loss features without the loss on the decoy site
        plain = df_candidates[is_decoy].copy()
        plain["SpecId"] = plain["SpecId"] + PLAIN_DECOY_SUFFIX
        for feature, plain_feature in PLAIN_DECOY_FEATURES.items():
            if feature in plain.columns:
                plain[feature] = candidates.loc[plain.index, plain_feature]
        df_candidates = pd.concat([df_candidates, plain], ignore_index=True)
        is_decoy = pd.concat([is_decoy, pd.Series(True, index=plain.index)], ignore_index=True)
    # deterministic order (percolator assigns PSMs to cross validation bins by input order)
    order = df_candidates["spectral_angle"].sort_values(ascending=False, kind="mergesort").index
    df_candidates, is_decoy = df_candidates.loc[order], is_decoy.loc[order]

    train_dir.mkdir(parents=True, exist_ok=True)
    df_candidates[~is_decoy].to_csv(train_dir / "rescore.tab", sep="\t", index=False)
    rescore_with_percolator(input_file=train_dir / "rescore.tab", output_folder=train_dir, num_threads=num_threads)

    score_dir.mkdir(parents=True, exist_ok=True)
    averaged_weights = score_dir / "averaged.weights.csv"
    _average_fold_weights(train_dir / "rescore.percolator.weights.csv", averaged_weights)
    constrained_weights = score_dir / "constrained.weights.csv"
    _constrain_weights(
        averaged_weights, constrained_weights, LOCALIZATION_POSITIVE_FEATURES + list(extra_positive_features or [])
    )
    df_candidates["ExpMass"] = np.arange(1, len(df_candidates) + 1)
    df_candidates.to_csv(score_dir / "rescore.tab", sep="\t", index=False)
    rescore_with_percolator(
        input_file=score_dir / "rescore.tab",
        output_folder=score_dir,
        init_weights=constrained_weights,
        num_threads=num_threads,
    )

    psms = pd.concat(
        [
            pd.read_csv(score_dir / "rescore.percolator.psms.txt", sep="\t", dtype={"filename": str}),
            pd.read_csv(score_dir / "rescore.percolator.decoy.psms.txt", sep="\t", dtype={"filename": str}),
        ],
        ignore_index=True,
    )
    psms["ScanNr"] = _scan_numbers(psms)
    spectrum = ["filename", "ScanNr"]
    # probabilities relative to the best candidate of the spectrum (incl. decoys)
    if probability == "score":
        # log odds of the PEP fit of the training step
        log_odds = _score_slope(train_dir / "rescore.percolator.psms.txt") * psms["score"]
    else:
        log_odds = -np.log(psms["posterior_error_prob"].clip(lower=1e-15))
    psms["score_prob"] = np.exp(log_odds - log_odds.groupby([psms[column] for column in spectrum]).transform("max"))
    decoy_columns = ["decoy_win", "decoy_weight", "decoy_win_plain", "decoy_weight_plain"]
    if decoy_pattern:
        decoys = psms["peptide"].str.contains(decoy_pattern, regex=True)
        plain = psms["PSMId"].str.endswith(PLAIN_DECOY_SUFFIX)
        real = psms[~decoys]
        best_real = real.groupby(spectrum)["score_prob"].max()
        flags = {}
        for suffix, decoy_type in (("", decoys & ~plain), ("_plain", plain if use_plain else decoys & ~plain)):
            typed = psms[~decoys | decoy_type]
            weight = _decoy_weights(typed, decoy_type[typed.index], ptm_pattern, decoy_weighting)
            best_decoy = psms[decoy_type].groupby(spectrum)["score_prob"].max()
            win = best_decoy.reindex(best_real.index).fillna(0) > best_real
            flags[f"decoy_win{suffix}"] = win
            flags[f"decoy_weight{suffix}"] = weight.reindex(best_real.index).fillna(0).where(win, 0.0)
        psms = real.join(pd.DataFrame(flags), on=spectrum)
    else:
        psms = psms.assign(decoy_win=False, decoy_weight=0.0, decoy_win_plain=False, decoy_weight_plain=0.0)
    psms["loc_prob"] = psms["score_prob"] / psms.groupby(spectrum)["score_prob"].transform("sum")
    psms = psms.sort_values(spectrum + ["loc_prob"], ascending=[True, True, False])
    loc_probs = psms.groupby(spectrum)["loc_prob"]
    second_best = loc_probs.transform(lambda probs: probs.iloc[1] if len(probs) > 1 else 0.0)
    psms["loc_prob_delta"] = loc_probs.transform("first") - second_best
    psms["possible_sites"] = loc_probs.transform("size")
    # model based error probabilities of every candidate: 1 - softmax of the scores of the real candidates of the
    # spectrum, sharper than the training slope (see MODEL_ERROR_SHARPNESS), used by the "hybrid" FLR estimator
    slope = _score_slope(train_dir / "rescore.percolator.psms.txt")
    keys = [psms[column] for column in spectrum]
    for column, sharpness in MODEL_ERROR_SHARPNESS.items():
        logits = sharpness * slope * psms["score"]
        maximum = logits.groupby(keys).transform("max")
        log_sum = np.log(np.exp(logits - maximum).groupby(keys).transform("sum")) + maximum
        psms[column] = 1 - np.exp(logits - log_sum)
    return psms[
        ["PSMId", "filename", "ScanNr"]
        + (["peptide", "score"] if keep_scores else [])
        + ["loc_prob", "loc_prob_delta", "possible_sites"]
        + decoy_columns
        + list(MODEL_ERROR_SHARPNESS)
    ]


def _best_localization(df_rescore: pd.DataFrame, localization: pd.DataFrame) -> pd.DataFrame:
    """Keep the localization candidate with the highest localization probability per spectrum."""
    best = localization.sort_values("loc_prob", ascending=False).drop_duplicates(["filename", "ScanNr"])
    return df_rescore[df_rescore["SpecId"].isin(best["PSMId"])]


def _sequence_and_charge(psm_ids: pd.Series, filenames: pd.Series, scans: pd.Series) -> tuple[list[str], list[int]]:
    """Modified sequence and precursor charge of PSM ids <filename>-<scan>-<sequence>-<charge>-<scan event>."""
    rest = [
        psm_id[len(f"{filename}-{scan}-") :] for psm_id, filename, scan in zip(psm_ids, filenames, scans, strict=True)
    ]
    parts = [r.rsplit("-", 2) for r in rest]
    return [p[0] for p in parts], [int(p[1]) for p in parts]


def _mod_prefix_counts(sequence: str) -> np.ndarray:
    """Number of every modification on residues 1..k of the sequence (rows k = 1..n, one column per unimod id)."""
    core = re.sub(r"^(?:\[UNIMOD:\d+\])+-", "", sequence)
    tokens = _TOKEN.findall(core)
    per_residue = [re.findall(r"UNIMOD:(\d+)", mods) for _, mods in tokens]
    ids = sorted({u for mods in per_residue for u in mods})
    counts = np.array([[mods.count(u) for u in ids] for mods in per_residue], dtype=int).reshape(len(tokens), len(ids))
    return np.cumsum(counts, axis=0)


def _differing_ions(sequence_a: str, sequence_b: str, charge: int) -> np.ndarray:
    """
    Site-determining ions of two localizations of the same peptide, in the Prosit order (y1+1, y1+2, y1+3, b1+1, ...).

    A b or y ion differs if it carries a different set of modifications in the two localizations; only fragment
    charges up to min(3, precursor charge) are considered.
    """
    prefix_a, prefix_b = _mod_prefix_counts(sequence_a), _mod_prefix_counts(sequence_b)
    n = len(prefix_a)
    if prefix_a.shape != prefix_b.shape or n < 2:
        return np.zeros(max(n - 1, 0) * 6, dtype=bool)
    differs = (prefix_a != prefix_b).any(axis=1)  # b_k differs for k = 1..n
    k = np.arange(1, n)
    b_differs, y_differs = differs[k - 1], differs[n - k - 1]  # y_k carries residues n-k+1..n: differs like b_(n-k)
    return np.repeat(np.concatenate([y_differs[:, None], b_differs[:, None]], axis=1), 3, axis=1).ravel() & (
        _charge_mask(n - 1, charge)
    )


def _charge_mask(n_positions: int, charge: int) -> np.ndarray:
    """Fragment charges up to min(3, precursor charge), in the Prosit order (y1+1, y1+2, y1+3, b1+1, ...)."""
    return np.tile(np.arange(1, 4) <= min(3, charge), 2 * n_positions)


def _stored_annotation(data_dir: Path, filename: str):
    """Predicted and annotated intensities of the stored library of a raw file, with the row of every PSM."""
    from ..data import Spectra

    library = Spectra.from_hdf5(data_dir / f"{filename}.mzml.pred.hdf5")
    obs = library.obs.reset_index(drop=True)
    rows = {
        (int(scan), sequence, int(charge)): i
        for i, (scan, sequence, charge) in enumerate(
            zip(obs["SCAN_NUMBER"], obs["MODIFIED_SEQUENCE"], obs["PRECURSOR_CHARGE"], strict=True)
        )
    }
    dense = [library.layers[name] for name in ("pred_int", "raw_int")]
    dense = [layer.toarray() if hasattr(layer, "toarray") else np.asarray(layer) for layer in dense]
    return rows, np.clip(dense[0], 0, None), dense[1]


def sdi_features(
    candidates: pd.DataFrame, first_pass: pd.DataFrame, data_dir: Path, decoy_pattern: Optional[str] = None
) -> pd.DataFrame:
    """
    Site-determining ion (SDI) features of all localization candidates (SDI_FEATURES).

    A and B are the two best real candidates of a spectrum by the first pass scores (localization_probabilities
    without these features); the competitor of A is B, of every other candidate (decoys included) A. The SDIs of a
    candidate c are the possible b / y ions (see _differing_ions) that differ between c and its competitor. From the
    stored predictions (pred_int) and annotation (raw_int, i.e. the matches at the configured fragment tolerance):
    sdi_pred_obs_frac = predicted intensity of the observed SDIs / of all SDIs of c, sdi_pred_missing = 1 - predicted
    intensity of the unobserved SDIs / of all ions of c, sdi_none = 1 if c has no predicted SDI (or no competitor or
    no stored prediction), then both other features are 1 (no penalty).

    :param candidates: rescore features of all localization candidates (SpecId, filename, ScanNr, Peptide)
    :param first_pass: localization_probabilities of the candidates without the SDI features (PSMId, score)
    :param data_dir: directory with the stored libraries <filename>.mzml.pred.hdf5
    :param decoy_pattern: regex matching decoy candidates in the Peptide column, None if no decoys are used
    :return: SpecId and SDI_FEATURES of every candidate
    """
    cand = candidates[["SpecId", "filename", "ScanNr", "Peptide"]].copy()
    cand["seq"], cand["charge"] = _sequence_and_charge(cand["SpecId"], cand["filename"], cand["ScanNr"])
    cand = cand.merge(first_pass[["PSMId", "score"]].rename(columns={"PSMId": "SpecId"}), on="SpecId", how="left")
    cand["decoy"] = cand["Peptide"].str.contains(decoy_pattern, regex=True) if decoy_pattern else False
    result = []
    for filename, file_candidates in cand.groupby("filename", sort=False):
        rows, pred, raw = _stored_annotation(data_dir, filename)
        for scan, group in file_candidates.groupby("ScanNr", sort=False):
            reals = group[~group["decoy"] & group["score"].notna()].sort_values("score", ascending=False)
            best = list(reals["SpecId"].iloc[:2])
            info = {}
            for spec_id, seq, charge in zip(group["SpecId"], group["seq"], group["charge"], strict=True):
                i = rows.get((int(scan), seq, charge))
                info[spec_id] = None if i is None else (seq, charge, i)
            a = best[0] if best else None
            b = best[1] if len(best) > 1 else None
            for spec_id in group["SpecId"]:
                competitor = b if spec_id == a else a
                c, q = info[spec_id], info.get(competitor) if competitor else None
                if c is None or q is None:
                    result.append((spec_id, 1.0, 1.0, 1.0))
                    continue
                differs = _differing_ions(c[0], q[0], c[1])
                n = len(differs)
                p, o, o_q = pred[c[2], :n], raw[c[2], :n], raw[q[2], :n]
                possible = (o > 0) & _charge_mask(n // 6, c[1])
                sdi = possible & (o_q > 0) & differs
                total, sdi_pred = p[possible].sum(), p[sdi].sum()
                if sdi_pred <= 0 or total <= 0:
                    result.append((spec_id, 1.0, 1.0, 1.0))
                    continue
                seen = sdi & (o > OBSERVED_INTENSITY)
                result.append((spec_id, float(p[seen].sum() / sdi_pred), float(1 - p[sdi & ~seen].sum() / total), 0.0))
    return pd.DataFrame(result, columns=["SpecId"] + SDI_FEATURES)


def ambiguity_flags(localization: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    """
    Flag localizations without site-specific evidence against the runner-up (ambiguity_check).

    The call A of a spectrum is its candidate with the highest localization probability, the runner-up B the
    highest scored other real candidate. From the stored annotation (raw_int, i.e. the matches at the configured
    fragment tolerance) of the site-determining ions of A and B (_differing_ions): no_sdi = none of them is observed
    for A or B, favours_B = the observed SDIs of B have more intensity than those of A (ions possible for both). The
    call is ambiguous if no_sdi or favours_B. Spectra with a single real candidate or without stored annotation are
    not flagged.

    :param localization: localization_probabilities output (PSMId, filename, ScanNr, peptide, score, loc_prob)
    :param data_dir: directory with the stored libraries <filename>.mzml.pred.hdf5
    :return: localization with localization_ambiguous (bool, on the calls) and localization_runner_up (peptide of B
        on the calls)
    """
    spectrum = ["filename", "ScanNr"]
    localization = localization.copy()
    called = localization.sort_values("loc_prob", ascending=False).drop_duplicates(spectrum)
    others = localization[~localization["PSMId"].isin(called["PSMId"])]
    runner_up = others.sort_values("score", ascending=False).drop_duplicates(spectrum)
    pairs = called[spectrum + ["PSMId", "peptide"]].merge(
        runner_up[spectrum + ["PSMId", "peptide"]], on=spectrum, suffixes=("", "_B")
    )
    pairs["seq"], pairs["charge"] = _sequence_and_charge(pairs["PSMId"], pairs["filename"], pairs["ScanNr"])
    pairs["seq_B"] = _sequence_and_charge(pairs["PSMId_B"], pairs["filename"], pairs["ScanNr"])[0]
    flagged = np.zeros(len(pairs), dtype=bool)
    for filename, group in pairs.groupby("filename", sort=False):
        rows, _, raw = _stored_annotation(data_dir, filename)
        for i, scan, seq_a, seq_b, charge in zip(
            group.index, group["ScanNr"], group["seq"], group["seq_B"], group["charge"], strict=True
        ):  # pairs has a RangeIndex: i is the position
            row_a, row_b = rows.get((int(scan), seq_a, charge)), rows.get((int(scan), seq_b, charge))
            if row_a is None or row_b is None:
                continue
            differs = _differing_ions(seq_a, seq_b, charge)
            n = len(differs)
            ra, rb = raw[row_a, :n], raw[row_b, :n]
            observed_a, observed_b = differs & (ra > OBSERVED_INTENSITY), differs & (rb > OBSERVED_INTENSITY)
            both = differs & (ra > 0) & (rb > 0)
            favours_b = rb[both & observed_b].sum() > ra[both & observed_a].sum()
            flagged[i] = not (observed_a.any() or observed_b.any()) or favours_b
    pairs["localization_ambiguous"] = flagged
    localization = localization.merge(
        pairs[["PSMId", "localization_ambiguous", "peptide_B"]].rename(columns={"peptide_B": "localization_runner_up"}),
        on="PSMId",
        how="left",
    )
    localization["localization_ambiguous"] = localization["localization_ambiguous"].eq(True)
    return localization


def localization_groups(peptides: pd.Series, ptm_pattern: str) -> pd.Series:
    r"""
    Get the residue of the called PTM sites of every PSM: S, T or Y, or "mixed" if the sites are on different residues.

    :param peptides: localized peptides in percolator format
    :param ptm_pattern: regex matching the localized PTMs, e.g. \[UNIMOD:(?:21)\]
    :return: the residue group of every PSM, empty for PSMs without the PTMs
    """
    residues = peptides.str.findall(rf"([A-Z]){ptm_pattern}").map(set)
    return residues.map(lambda found: next(iter(found)) if len(found) == 1 else ("mixed" if found else ""))


def neighbor_classes(peptides: pd.Series, ptm_pattern: str, possible_residues: str) -> pd.Series:
    """
    Get whether a possible site that is not localized is directly next to a localized site of every PSM.

    Wrong localizations are mostly on a neighboring site, so these localizations have a higher FLR at the same score.

    :param peptides: localized peptides in percolator format
    :param ptm_pattern: regex matching the localized PTMs
    :param possible_residues: residues the PTMs can be located on, e.g. "STY"
    :return: "neighbor" or "no neighbor" for every PSM
    """
    core = peptides.str.replace(r"^[^.]*\.|\.[^.]*$", "", regex=True).str.replace(
        r"^(?:\[UNIMOD:\d+\])+-", "", regex=True
    )
    sites = _site_positions(core, ptm_pattern)
    residues = core.str.replace(r"\[UNIMOD:\d+\]", "", regex=True)

    def has_neighbor(sequence: str, positions: frozenset) -> bool:
        return any(
            0 < q <= len(sequence) and q not in positions and sequence[q - 1] in possible_residues
            for p in positions
            for q in (p - 1, p + 1)
        )

    return pd.Series(
        ["neighbor" if has_neighbor(seq, pos) else "no neighbor" for seq, pos in zip(residues, sites, strict=True)],
        index=peptides.index,
    )


def _ranked_flr(psms: pd.DataFrame) -> pd.Series:
    """
    FLR q-values of PSMs ranked by loc_prob_delta, then loc_prob: cumulative decoy_weight / number of PSMs.

    With the column localization_ambiguous (ambiguity_check), ambiguous localizations are ranked after all others.
    """
    if "localization_ambiguous" in psms.columns:
        ranked = psms.assign(_ambiguous=psms["localization_ambiguous"].eq(True)).sort_values(
            ["_ambiguous", "loc_prob_delta", "loc_prob"], ascending=[True, False, False], kind="mergesort"
        )
    else:
        ranked = psms.sort_values(["loc_prob_delta", "loc_prob"], ascending=False, kind="mergesort")
    flr = ranked["decoy_weight"].fillna(0).cumsum().values / np.arange(1, len(ranked) + 1)
    return pd.Series(np.minimum.accumulate(flr[::-1])[::-1], index=ranked.index)


def _model_flr(psms: pd.DataFrame, error_column: str, group_column: str) -> pd.Series:
    """
    FLR q-values from model based error probabilities, per group.

    Within every group, PSMs are ranked by their error probability and the FLR of the top k PSMs is the mean error
    probability of these PSMs, reported as q-value. With the column localization_ambiguous (ambiguity_check),
    ambiguous localizations are ranked after all others of their group, each part by its error probability.
    """
    flr = pd.Series(np.nan, index=psms.index)
    ambiguous = "localization_ambiguous" in psms.columns
    for _, group in psms.groupby(group_column):
        if ambiguous:
            order = group.assign(_ambiguous=group["localization_ambiguous"].eq(True)).sort_values(
                ["_ambiguous", error_column], kind="mergesort"
            )
            ranked = order[error_column]
        else:
            ranked = group[error_column].sort_values(kind="mergesort")
        running = ranked.cumsum().to_numpy() / np.arange(1, len(ranked) + 1)
        flr.loc[ranked.index] = np.minimum.accumulate(running[::-1])[::-1]
    return flr


def hybrid_localization_flr(psms: pd.DataFrame, group_column: str, fdr: float = 0.01) -> pd.DataFrame:
    """
    Combine the decoy estimated FLR with model based FLRs, so that the estimate is reliable above 1-2% FLR.

    Decoy candidates only reveal wrong localizations if a decoy scores above the called site, which rarely happens
    once a wrong site has already beaten the true one, so the decoy estimate saturates for localizations with a
    possible site next to the called site. The model based FLR uses the error probability of the call, 1 - softmax of
    the scores of the real candidates of the spectrum, summed within groups of the called residue and whether a
    possible site is next to the called site. localization_flr = max(model FLR (MODEL_SHARPNESS_HIGH),
    min(decoy FLR, model FLR (MODEL_SHARPNESS_LOW))): at the low end, the less conservative of the decoy and the
    model estimate, above about 1-2% the more conservative one. The decoy estimate is kept as localization_flr_decoy.

    With the column localization_ambiguous (ambiguity_check), the ambiguous localizations are ranked after all others
    of their group, and their FLR is at least the decoy FLR estimated over the ambiguous localizations alone (pooled
    over all groups): the model errors of calls without site-specific evidence are too optimistic, while their decoys
    still compete on equal terms.

    :param psms: PSMs with localization_flr (decoy estimate), model_error_high, model_error_low and group_column
    :param group_column: column with the group of the model based FLR
    :param fdr: PSM level q-value cutoff for the PSMs used to estimate the FLR
    :return: psms with localization_flr replaced by the hybrid estimate and localization_flr_decoy
    """
    psms["localization_flr_decoy"] = psms["localization_flr"]
    passing = psms[(psms["q-value"] <= fdr) & psms["loc_prob"].notna()]
    if passing.empty:
        return psms
    high = _model_flr(passing, "model_error_high", group_column)
    low = _model_flr(passing, "model_error_low", group_column)
    decoy = passing["localization_flr"]
    hybrid = np.maximum(high, np.minimum(decoy, low))
    if "localization_ambiguous" in passing.columns:
        ambiguous = passing[passing["localization_ambiguous"].eq(True)]
        if not ambiguous.empty:
            hybrid.loc[ambiguous.index] = np.maximum(hybrid.loc[ambiguous.index], _ranked_flr(ambiguous))
    psms.loc[passing.index, "localization_flr"] = hybrid
    return psms


def estimate_localization_flr(
    psms: pd.DataFrame, fdr: float = 0.01, group_column: Optional[str] = None, min_group_size: int = 200
) -> pd.DataFrame:
    """
    Estimate the false localization rate (FLR) from decoy residues.

    PSMs (q-value <= fdr) are ranked by loc_prob_delta, then loc_prob. In every spectrum flagged as decoy_win, a
    candidate with the PTM on a decoy residue (e.g. phosphorylation on A) scored higher than all real candidates. Each
    decoy win stands for decoy_weight expected false localizations among the real candidates of that spectrum, with
    decoy and wrong real candidates matched by their distance to the called site (see _decoy_weights). FLR = sum of
    decoy_weight / number of localizations, reported as q-value (localization_flr) for every PSM up to its rank.

    With group_column (e.g. the residue of the called sites, see localization_groups), the FLR is estimated within
    every group of at least min_group_size PSMs, since the decoys can represent the wrong localizations of the groups
    differently well (e.g. phosphorylated S/T lose phosphoric acid, Y does not). Smaller groups keep the FLR estimated
    over all PSMs.

    :param psms: percolator PSMs with loc_prob, loc_prob_delta, decoy_win and decoy_weight
    :param fdr: PSM level q-value cutoff for the PSMs used to estimate the FLR
    :param group_column: column with the group of every PSM, None to estimate the FLR over all PSMs
    :param min_group_size: smallest number of PSMs of a group with its own FLR estimate
    :return: psms with the column localization_flr (NaN for PSMs above the fdr)
    """
    psms["localization_flr"] = np.nan
    passing = psms[(psms["q-value"] <= fdr) & psms["loc_prob"].notna()]
    if passing.empty:
        return psms
    psms.loc[passing.index, "localization_flr"] = _ranked_flr(passing)
    if group_column is not None:
        for _, group in passing.groupby(group_column):
            if len(group) >= min_group_size:
                psms.loc[group.index, "localization_flr"] = _ranked_flr(group)
    return psms


def _calibrate(reference: pd.DataFrame, deltas: pd.Series, anchors: tuple) -> np.ndarray:
    """
    Interpolate the localization score of loc_prob_delta values from the FLR anchors.

    :param reference: PSMs with localization_flr and loc_prob_delta that define the delta cutoffs of the anchors
    :param deltas: the loc_prob_delta values to calibrate
    :param anchors: pairs of (estimated FLR, localization score)
    :return: localization scores
    """
    cutoffs, scores = [1.0], [1.0]
    for flr, score in sorted(anchors):
        within = reference.loc[reference["localization_flr"] <= flr, "loc_prob_delta"]
        # the delta cutoff of this FLR; nothing reaches the score if no PSM is within the FLR
        cutoff = float(within.min()) if not within.empty else 1.0
        cutoffs.append(min(cutoff, cutoffs[-1]))
        scores.append(score)
    cutoffs.append(0.0)
    scores.append(0.0)
    cutoffs, scores = np.array(cutoffs[::-1]), np.array(scores[::-1])
    # strictly increasing cutoffs for the interpolation
    cutoffs = cutoffs + np.arange(len(cutoffs)) * 1e-12
    return np.interp(deltas, cutoffs, scores)


def calibrate_localization_score(
    psms: pd.DataFrame, anchors: tuple = ((0.01, 0.75), (0.05, 0.5)), group_column: Optional[str] = None
) -> pd.DataFrame:
    """
    Calibrate loc_prob_delta to a localization score using the decoy estimated FLR.

    For every anchor (FLR, score), the loc_prob_delta at which the estimated localization_flr reaches FLR is mapped
    to score, e.g. a localization score of 0.75 corresponds to 1% and 0.5 to 5% estimated FLR. The score is linearly
    interpolated between the anchors, 0 and 1. With group_column, every group is calibrated separately (e.g. per
    residue of the called sites, with the FLR estimated per group).

    :param psms: percolator PSMs with loc_prob_delta and localization_flr
    :param anchors: pairs of (estimated FLR, localization score), with decreasing score for increasing FLR
    :param group_column: column with the group of every PSM, None to calibrate all PSMs together
    :return: psms with the column localization_score (NaN for PSMs without loc_prob_delta)
    """
    psms["localization_score"] = np.nan
    passing = psms[psms["localization_flr"].notna()]
    if passing.empty:
        return psms
    has_delta = psms["loc_prob_delta"].notna()
    if group_column is None:
        psms.loc[has_delta, "localization_score"] = _calibrate(passing, psms.loc[has_delta, "loc_prob_delta"], anchors)
        return psms
    for group, members in psms[has_delta].groupby(group_column):
        reference = passing[passing[group_column] == group]
        psms.loc[members.index, "localization_score"] = _calibrate(reference, members["loc_prob_delta"], anchors)
    return psms


def _add_localization_probabilities(
    result_dir: Path,
    localization: pd.DataFrame,
    use_decoys: bool,
    ptm_pattern: str = r"\[UNIMOD:21\]",
    possible_residues: Optional[str] = None,
    flr_groups: str = "none",
    flr_estimator: str = "hybrid",
    ptm_residues: Optional[str] = None,
):
    """
    Write the percolator results of result_dir with localization probabilities (*_with_prob.*psms.txt).

    If decoy residues are used, the target PSMs also get the estimated localization_flr and the calibrated
    localization_score. With flr_groups "none" (default) they are estimated over all PSMs, with "residue" per residue
    of the called sites. The called residues in PLAIN_DECOY_GROUPS use the plain decoys (decoy_win_plain,
    decoy_weight_plain) if the decoys are annotated with the neutral loss, the others the decoys as annotated. With
    possible_residues, the groups are further split by whether a possible site is next to the called site
    (neighbor_classes). With flr_estimator "hybrid", the decoy estimate is combined with model based estimates per
    called residue and neighbor class of ptm_residues (see hybrid_localization_flr), with the sharper model errors for
    calls without a neighboring possible site (see MODEL_ERROR_SHARPNESS). If the localization has the column
    localization_ambiguous (ambiguity_check), the ambiguous localizations are ranked after all others of their group
    by both estimates, and their FLR is at least the decoy FLR of the ambiguous localizations alone.
    """
    columns = ["PSMId", "loc_prob", "loc_prob_delta", "possible_sites", "decoy_win", "decoy_weight"]
    columns += ["decoy_win_plain", "decoy_weight_plain"] + list(MODEL_ERROR_SHARPNESS)
    # ambiguity_check: ambiguous localizations are ranked after all others by both FLR estimates (see ambiguity_flags)
    columns += [column for column in ("localization_ambiguous", "localization_runner_up") if column in localization]
    for psm_type in ["psms", "decoy.psms"]:
        psms = pd.read_csv(result_dir / f"rescore.percolator.{psm_type}.txt", sep="\t", dtype={"filename": str})
        psms = psms.merge(localization[columns], on="PSMId", how="left")
        if "localization_ambiguous" in psms.columns:
            psms["localization_ambiguous"] = psms["localization_ambiguous"].eq(True)
        if use_decoys and psm_type == "psms":
            residues = localization_groups(psms["peptide"], ptm_pattern)
            plain = residues.isin(PLAIN_DECOY_GROUPS)
            psms["decoy_win"] = psms["decoy_win"].where(~plain, psms["decoy_win_plain"])
            psms["decoy_weight"] = psms["decoy_weight"].where(~plain, psms["decoy_weight_plain"])
            groups = [residues] if flr_groups == "residue" else []
            if possible_residues:
                groups.append(neighbor_classes(psms["peptide"], ptm_pattern, possible_residues))
            group_column = None
            if groups:
                psms["localization_group"] = groups[0] if len(groups) == 1 else groups[0] + " " + groups[1]
                group_column = "localization_group"
            psms = estimate_localization_flr(psms, group_column=group_column)
            # the model based estimate needs the residues the PTMs can be located on for the neighbor classes
            if flr_estimator == "hybrid" and ptm_residues:
                neighbors = neighbor_classes(psms["peptide"], ptm_pattern, ptm_residues)
                psms["model_group"] = residues + " " + neighbors
                for column in ("model_error_high", "model_error_low"):
                    psms[column] = psms[column].where(neighbors == "neighbor", psms[f"{column}_no_neighbor"])
                psms = hybrid_localization_flr(psms, group_column="model_group")
            psms = calibrate_localization_score(psms, group_column=group_column)
        psms = psms.drop(columns=["model_error_high_no_neighbor", "model_error_low_no_neighbor"])
        psms.to_csv(result_dir / f"rescore.percolator_with_prob.{psm_type}.txt", sep="\t", index=False)


def _localize_candidates(
    candidates: pd.DataFrame,
    fdr_dir: Path,
    decoy_pattern: Optional[str],
    ptm_pattern: str,
    probability: str,
    decoy_weighting: str,
    data_dir: Optional[Path],
    use_sdi_features: bool,
    ambiguity_check: bool,
    num_threads: int = 3,
) -> pd.DataFrame:
    """
    Score the localization candidates (localization_probabilities in localize_mod / localize_scores).

    With use_sdi_features, the candidates are first scored without the SDI features (localize_mod_first_pass /
    localize_scores_first_pass) to find the two best candidates of every spectrum, then scored with SDI_FEATURES
    (see sdi_features). With ambiguity_check, the calls are flagged by ambiguity_flags. Both need the stored libraries
    in data_dir. Percolator runs with num_threads threads.
    """
    if (use_sdi_features or ambiguity_check) and data_dir is None:
        raise ValueError("sdi_features and ambiguity_check need the directory of the stored libraries (data_dir).")
    extra, extra_positive = None, None
    if use_sdi_features:
        first_pass = localization_probabilities(
            candidates,
            fdr_dir / "localize_mod_first_pass",
            fdr_dir / "localize_scores_first_pass",
            decoy_pattern,
            ptm_pattern,
            probability,
            decoy_weighting,
            keep_scores=True,
            num_threads=num_threads,
        )
        features = sdi_features(candidates, first_pass, data_dir, decoy_pattern)
        candidates = candidates.merge(features, on="SpecId", how="left")
        candidates[SDI_FEATURES] = candidates[SDI_FEATURES].fillna(1.0)
        extra, extra_positive = SDI_FEATURES, SDI_POSITIVE_FEATURES
    localization = localization_probabilities(
        candidates,
        fdr_dir / "localize_mod",
        fdr_dir / "localize_scores",
        decoy_pattern,
        ptm_pattern,
        probability,
        decoy_weighting,
        extra,
        extra_positive,
        keep_scores=ambiguity_check,
        num_threads=num_threads,
    )
    if ambiguity_check:
        localization = ambiguity_flags(localization, data_dir)
    return localization


def localize_ptms(
    fdr_dir: Path,
    ptm_sites: dict[int, list[str]],
    decoy_residues: Optional[list[str]] = None,
    neighbor_split: bool = False,
    probability: str = "pep",
    decoy_weighting: str = "count",
    flr_groups: str = "none",
    flr_estimator: str = "hybrid",
    use_sdi_features: bool = False,
    ambiguity_check: bool = False,
    data_dir: Optional[Path] = None,
    num_threads: int = 3,
):
    """
    Rescore with percolator and localize PTMs.

    Citrullination (UNIMOD:7) and phosphorylation (UNIMOD:21) have dedicated workflows if they are the only PTM to
    localize. Any other PTM or combination of PTMs is localized with the generic workflow.

    :param fdr_dir: the output directory containing rescore.tab with all localization candidates
    :param ptm_sites: the localized PTMs {unimod id: [residues]}
    :param decoy_residues: decoy residues used to estimate the false localization rate, None if not used
    :param neighbor_split: estimate the FLR separately for localizations with and without a possible site next to
        the called site
    :param probability: localization probabilities from the PEPs ("pep") or the scores ("score") of all candidates,
        see localization_probabilities
    :param decoy_weighting: expected false localizations of a decoy win, "count" or "distance", see _decoy_weights
    :param flr_groups: estimate the FLR over all PSMs ("none") or per residue of the called sites ("residue")
    :param flr_estimator: FLR from the decoys only ("decoy") or combined with model based estimates ("hybrid"), see
        hybrid_localization_flr
    :param use_sdi_features: score the candidates with site-determining ion features (SDI_FEATURES, see sdi_features)
    :param ambiguity_check: rank localizations without site-specific evidence against the runner-up after all others
        when estimating the FLR (see ambiguity_flags)
    :param data_dir: directory with the stored libraries <raw file>.mzml.pred.hdf5, needed by use_sdi_features and
        ambiguity_check
    :param num_threads: number of threads used by percolator
    """
    df_rescore = pd.read_csv(fdr_dir / "rescore.tab", sep="\t", dtype={"filename": str})
    unimod_id = next(iter(ptm_sites)) if len(ptm_sites) == 1 else 0
    ptm_pattern = rf"UNIMOD:(?:{_ptm_ids(ptm_sites)})\]"
    decoy_pattern = _decoy_site_pattern(ptm_sites, decoy_residues)
    use_decoys = decoy_pattern is not None
    ptm_residues = "".join(sorted({aa for sites in ptm_sites.values() for aa in sites}))
    possible_residues = ptm_residues if neighbor_split else None

    new_rescore_dir = fdr_dir / "localize_mod"
    new_rescore_dir.mkdir(parents=True, exist_ok=True)

    if unimod_id == 7:
        df_rescore["id"] = df_rescore["filename"] + df_rescore["ScanNr"].astype(str)
        df_rescore_fil = df_rescore[df_rescore["Peptide"].str.contains("UNIMOD:" + str(unimod_id))]
        df_rescore = df_rescore[df_rescore["id"].isin(df_rescore_fil["id"].tolist())]
        df_rescore.drop(columns=["id"], inplace=True)
        df_rescore.to_csv(fdr_dir / "rescore.tab", sep="\t", index=False)
        rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir, num_threads=num_threads)
        df_rescore_psms_targets = pd.read_csv(
            fdr_dir / "rescore.percolator.psms.txt", sep="\t", dtype={"filename": str}
        )
        df_rescore_psms_decoys = pd.read_csv(
            fdr_dir / "rescore.percolator.decoy.psms.txt", sep="\t", dtype={"filename": str}
        )
        df_rescore_psms_targets = df_rescore_psms_targets[
            df_rescore_psms_targets["peptide"].apply(lambda x: "UNIMOD:" + str(unimod_id) in x)
        ]
        df_rescore_psms_decoys = df_rescore_psms_decoys[
            df_rescore_psms_decoys["peptide"].apply(lambda x: "UNIMOD:" + str(unimod_id) in x)
        ]
        df_rescore_psms = pd.concat([df_rescore_psms_targets, df_rescore_psms_decoys])
        df_rescore_psms = df_rescore_psms[["PSMId"]]
        df_rescore_psms.rename(columns={"PSMId": "SpecId"}, inplace=True)
        df_rescore = df_rescore.merge(df_rescore_psms, on="SpecId", how="inner")

        remove_cols = []
        for col in df_rescore.columns:
            if "_mod" in col:
                remove_cols.append(col)
        df_rescore.drop(columns=remove_cols, inplace=True)

        df_rescore.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
        rescore_with_percolator(
            input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir, num_threads=num_threads
        )
    elif unimod_id == 21:
        df_rescore = df_rescore[df_rescore["Peptide"].str.contains("UNIMOD:21")]

        # score all localization candidates and keep the most likely localization per spectrum
        localization = _localize_candidates(
            df_rescore,
            fdr_dir,
            decoy_pattern,
            rf"\[{ptm_pattern}",
            probability,
            decoy_weighting,
            data_dir,
            use_sdi_features,
            ambiguity_check,
            num_threads,
        )
        df_rescore = _best_localization(df_rescore, localization)

        if "P" in ptm_sites[unimod_id] or "H" in ptm_sites[unimod_id]:
            for site in ptm_sites[unimod_id]:
                new_rescore_dir = fdr_dir / ("localize_rt_all_" + site)
                new_rescore_dir.mkdir(parents=True, exist_ok=True)
                df_rescore_fil = df_rescore[df_rescore["Peptide"].str.contains(site + r"\[")]
                df_rescore_fil.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
                rescore_with_percolator(
                    input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir, num_threads=num_threads
                )
                _add_localization_probabilities(
                    new_rescore_dir,
                    localization,
                    use_decoys,
                    rf"\[{ptm_pattern}",
                    possible_residues,
                    flr_groups,
                    flr_estimator,
                    ptm_residues,
                )

        else:
            new_rescore_dir = fdr_dir / ("localize_phospho")
            new_rescore_dir.mkdir(parents=True, exist_ok=True)

            df_rescore.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
            rescore_with_percolator(
                input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir, num_threads=num_threads
            )
            _add_localization_probabilities(
                new_rescore_dir,
                localization,
                use_decoys,
                rf"\[{ptm_pattern}",
                possible_residues,
                flr_groups,
                flr_estimator,
                ptm_residues,
            )

    else:
        df_rescore_all = df_rescore
        df_rescore = df_rescore[df_rescore["Peptide"].str.contains(ptm_pattern, regex=True)]
        # score all localization candidates and keep the most likely localization per spectrum
        localization = _localize_candidates(
            df_rescore,
            fdr_dir,
            decoy_pattern,
            rf"\[{ptm_pattern}",
            probability,
            decoy_weighting,
            data_dir,
            use_sdi_features,
            ambiguity_check,
        )
        df_rescore_psms = _best_localization(df_rescore, localization)[["SpecId"]]
        # keep the localized PSMs of the PTMs and all PSMs without them
        df_rescore = df_rescore_all[
            df_rescore_all["SpecId"].isin(df_rescore_psms["SpecId"].values)
            | ~df_rescore_all["Peptide"].str.contains(ptm_pattern, regex=True)
        ]
        new_rescore_dir = fdr_dir / "mod_fdr"
        new_rescore_dir.mkdir(parents=True, exist_ok=True)

        df_rescore.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
        rescore_with_percolator(
            input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir, num_threads=num_threads
        )
        _add_localization_probabilities(
            new_rescore_dir,
            localization,
            use_decoys,
            rf"\[{ptm_pattern}",
            possible_residues,
            flr_groups,
            flr_estimator,
            ptm_residues,
        )
