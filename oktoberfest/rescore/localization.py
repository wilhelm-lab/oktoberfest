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
    "calibrate_localization_score",
    "localization_groups",
    "neighbor_classes",
    "decoy_proxy_fractions",
    "add_adjacent_decoy_candidates",
    "localize_ptms",
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
    # site determining ions: the ions that distinguish a localization from the closest alternative site
    "fraction_observed_sdi",
    "fraction_observed_and_predicted_sdi",
    "fraction_not_observed_but_predicted_sdi",
    "observed_intensity_sdi",
    "sum_observed_and_predicted_sdi",
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
    "fraction_observed_sdi",
    "fraction_observed_and_predicted_sdi",
    "observed_intensity_sdi",
    "sum_observed_and_predicted_sdi",
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


def _constrain_weights(weights_file: Path, output_file: Path):
    """
    Set negative percolator weights of LOCALIZATION_POSITIVE_FEATURES to 0.

    The weights file contains, for every cross validation bin, the feature names, the normalized and the raw weights.

    :param weights_file: percolator weights file
    :param output_file: path of the constrained weights file
    """
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
                    str(max(0.0, float(value))) if name in LOCALIZATION_POSITIVE_FEATURES else value
                    for name, value in zip(names, fields, strict=True)
                )
            )
    output_file.write_text("\n".join(constrained) + "\n")


def _scan_numbers(psms: pd.DataFrame) -> list[int]:
    """Get scan numbers from the PSMId (<raw file>-<scan number>-...) of percolator results."""
    return [
        int(psm_id[len(raw_file) :].split("-")[1])
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


def _decoy_weights(psms: pd.DataFrame, decoys: pd.Series, ptm_pattern: str) -> pd.Series:
    """
    Expected number of false localizations per spectrum in which a decoy candidate wins, matched by site distance.

    Wrong localizations are most often next to the true site, decoy sites can be anywhere in the peptide. The chance
    q(d) that a candidate at distance d to the called sites scores above the called candidate is estimated from the
    decoy candidates of all spectra. A decoy win stands for sum(q(d) of the wrong real candidates) / sum(q(d) of the
    decoy candidates) expected false localizations in that spectrum.

    :param psms: all scored candidates (real and decoy) with filename, ScanNr, peptide and score_prob
    :param decoys: mask of the decoy candidates
    :param ptm_pattern: regex matching the localized PTMs
    :return: per spectrum (filename, ScanNr) the expected false localizations of a decoy win
    """
    spectrum = ["filename", "ScanNr"]
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
) -> pd.DataFrame:
    """
    Score all localization candidates of each spectrum and compute localization probabilities.

    First, percolator is trained on localization specific features of all candidates, where the candidates of a
    spectrum compete against each other (train_dir). Second, all candidates are scored with the trained weights
    (percolator --static), with a unique ExpMass per candidate so that percolator keeps all candidates of a spectrum
    instead of only the best one (score_dir). The localization probability of a candidate is its probability
    normalized over all candidates of the same spectrum, with the log odds a + b * score of the percolator PEP fit of
    the training step: loc_prob = exp(b * score) / sum(exp(b * score)), i.e. (1 / PEP) / sum(1 / PEP) for PEP << 1.
    The PEPs of the scoring step are not used, since percolator's PEP estimate on all candidates of all spectra can
    collapse to a step function for large inputs, making the candidates of a spectrum indistinguishable.

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
    :return: one row per (non decoy) candidate with PSMId, filename, ScanNr, loc_prob, loc_prob_delta (difference
        between the best and second best candidate of the spectrum), possible_sites (number of candidates of the
        spectrum), decoy_win (a decoy candidate scores higher than all candidates of the spectrum) and decoy_weight,
        and decoy_win_plain and decoy_weight_plain for the plain decoys (same as decoy_win / decoy_weight without them)
    """
    columns = (
        ["SpecId", "Label", "ScanNr", "filename", "ExpMass"]
        + [feature for feature in LOCALIZATION_FEATURES if feature in candidates.columns]
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
    rescore_with_percolator(input_file=train_dir / "rescore.tab", output_folder=train_dir)

    score_dir.mkdir(parents=True, exist_ok=True)
    averaged_weights = score_dir / "averaged.weights.csv"
    _average_fold_weights(train_dir / "rescore.percolator.weights.csv", averaged_weights)
    constrained_weights = score_dir / "constrained.weights.csv"
    _constrain_weights(averaged_weights, constrained_weights)
    df_candidates["ExpMass"] = np.arange(1, len(df_candidates) + 1)
    df_candidates.to_csv(score_dir / "rescore.tab", sep="\t", index=False)
    rescore_with_percolator(
        input_file=score_dir / "rescore.tab",
        output_folder=score_dir,
        init_weights=constrained_weights,
    )

    psms = pd.concat(
        [
            pd.read_csv(score_dir / "rescore.percolator.psms.txt", sep="\t"),
            pd.read_csv(score_dir / "rescore.percolator.decoy.psms.txt", sep="\t"),
        ],
        ignore_index=True,
    )
    psms["ScanNr"] = _scan_numbers(psms)
    spectrum = ["filename", "ScanNr"]
    # probabilities relative to the best candidate of the spectrum (incl. decoys), from the log odds of the training
    log_odds = _score_slope(train_dir / "rescore.percolator.psms.txt") * psms["score"]
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
            weight = _decoy_weights(typed, decoy_type[typed.index], ptm_pattern)
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
    return psms[["PSMId", "filename", "ScanNr", "loc_prob", "loc_prob_delta", "possible_sites"] + decoy_columns]


def _best_localization(df_rescore: pd.DataFrame, localization: pd.DataFrame) -> pd.DataFrame:
    """Keep the localization candidate with the highest localization probability per spectrum."""
    best = localization.sort_values("loc_prob", ascending=False).drop_duplicates(["filename", "ScanNr"])
    return df_rescore[df_rescore["SpecId"].isin(best["PSMId"])]


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
    """FLR q-values of PSMs ranked by loc_prob_delta, then loc_prob: cumulative decoy_weight / number of PSMs."""
    ranked = psms.sort_values(["loc_prob_delta", "loc_prob"], ascending=False, kind="mergesort")
    flr = ranked["decoy_weight"].fillna(0).cumsum().values / np.arange(1, len(ranked) + 1)
    return pd.Series(np.minimum.accumulate(flr[::-1])[::-1], index=ranked.index)


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
):
    """
    Write the percolator results of result_dir with localization probabilities (*_with_prob.*psms.txt).

    If decoy residues are used, the target PSMs also get the estimated localization_flr and the calibrated
    localization_score, both per residue of the called sites (localization_group). The groups in PLAIN_DECOY_GROUPS
    use the plain decoys (decoy_win_plain, decoy_weight_plain), the others the decoys as annotated. With
    possible_residues, the groups are further split by whether a possible site is next to the called site
    (neighbor_classes).
    """
    columns = ["PSMId", "loc_prob", "loc_prob_delta", "possible_sites", "decoy_win", "decoy_weight"]
    columns += ["decoy_win_plain", "decoy_weight_plain"]
    for psm_type in ["psms", "decoy.psms"]:
        psms = pd.read_csv(result_dir / f"rescore.percolator.{psm_type}.txt", sep="\t")
        psms = psms.merge(localization[columns], on="PSMId", how="left")
        if use_decoys and psm_type == "psms":
            psms["localization_group"] = localization_groups(psms["peptide"], ptm_pattern)
            plain = psms["localization_group"].isin(PLAIN_DECOY_GROUPS)
            psms["decoy_win"] = psms["decoy_win"].where(~plain, psms["decoy_win_plain"])
            psms["decoy_weight"] = psms["decoy_weight"].where(~plain, psms["decoy_weight_plain"])
            if possible_residues:
                classes = neighbor_classes(psms["peptide"], ptm_pattern, possible_residues)
                psms["localization_group"] = psms["localization_group"] + " " + classes
            psms = estimate_localization_flr(psms, group_column="localization_group")
            psms = calibrate_localization_score(psms, group_column="localization_group")
        psms.to_csv(result_dir / f"rescore.percolator_with_prob.{psm_type}.txt", sep="\t", index=False)


def localize_ptms(
    fdr_dir: Path,
    ptm_sites: dict[int, list[str]],
    decoy_residues: Optional[list[str]] = None,
    neighbor_split: bool = False,
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
    """
    df_rescore = pd.read_csv(fdr_dir / "rescore.tab", sep="\t")
    unimod_id = next(iter(ptm_sites)) if len(ptm_sites) == 1 else 0
    ptm_pattern = rf"UNIMOD:(?:{_ptm_ids(ptm_sites)})\]"
    decoy_pattern = _decoy_site_pattern(ptm_sites, decoy_residues)
    use_decoys = decoy_pattern is not None
    possible_residues = (
        "".join(sorted({aa for sites in ptm_sites.values() for aa in sites})) if neighbor_split else None
    )

    new_rescore_dir = fdr_dir / "localize_mod"
    new_rescore_dir.mkdir(parents=True, exist_ok=True)

    if unimod_id == 7:
        df_rescore["id"] = df_rescore["filename"] + df_rescore["ScanNr"].astype(str)
        df_rescore_fil = df_rescore[df_rescore["Peptide"].str.contains("UNIMOD:" + str(unimod_id))]
        df_rescore = df_rescore[df_rescore["id"].isin(df_rescore_fil["id"].tolist())]
        df_rescore.drop(columns=["id"], inplace=True)
        df_rescore.to_csv(fdr_dir / "rescore.tab", sep="\t", index=False)
        rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
        df_rescore_psms_targets = pd.read_csv(fdr_dir / "rescore.percolator.psms.txt", sep="\t")
        df_rescore_psms_decoys = pd.read_csv(fdr_dir / "rescore.percolator.decoy.psms.txt", sep="\t")
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
        rescore_with_percolator(input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir)
    elif unimod_id == 21:
        df_rescore = df_rescore[df_rescore["Peptide"].str.contains("UNIMOD:21")]

        # score all localization candidates and keep the most likely localization per spectrum
        localization = localization_probabilities(
            df_rescore, new_rescore_dir, fdr_dir / "localize_scores", decoy_pattern, rf"\[{ptm_pattern}"
        )
        df_rescore = _best_localization(df_rescore, localization)

        if "P" in ptm_sites[unimod_id] or "H" in ptm_sites[unimod_id]:
            for site in ptm_sites[unimod_id]:
                new_rescore_dir = fdr_dir / ("localize_rt_all_" + site)
                new_rescore_dir.mkdir(parents=True, exist_ok=True)
                df_rescore_fil = df_rescore[df_rescore["Peptide"].str.contains(site + r"\[")]
                df_rescore_fil.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
                rescore_with_percolator(input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir)
                _add_localization_probabilities(
                    new_rescore_dir, localization, use_decoys, rf"\[{ptm_pattern}", possible_residues
                )

        else:
            new_rescore_dir = fdr_dir / ("localize_phospho")
            new_rescore_dir.mkdir(parents=True, exist_ok=True)

            df_rescore.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
            rescore_with_percolator(input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir)
            _add_localization_probabilities(
                new_rescore_dir, localization, use_decoys, rf"\[{ptm_pattern}", possible_residues
            )

    else:
        df_rescore_all = df_rescore
        df_rescore = df_rescore[df_rescore["Peptide"].str.contains(ptm_pattern, regex=True)]
        # score all localization candidates and keep the most likely localization per spectrum
        localization = localization_probabilities(
            df_rescore, new_rescore_dir, fdr_dir / "localize_scores", decoy_pattern, rf"\[{ptm_pattern}"
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
        rescore_with_percolator(input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir)
        _add_localization_probabilities(
            new_rescore_dir, localization, use_decoys, rf"\[{ptm_pattern}", possible_residues
        )
