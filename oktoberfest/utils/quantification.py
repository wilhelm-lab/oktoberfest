import json
import logging
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable, Optional

import anndata
import numpy as np
import pandas as pd
from picked_group_fdr.digestion_params import DigestionParams
from picked_group_fdr.pipeline import pipeline
from pyteomics import mass, mzml

from . import Config

logger = logging.getLogger(__name__)

DECOY_PREFIX = "rev_"
SPEC_ID_COLUMNS = ["RAW_FILE", "SCAN_NUMBER", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"]
# ppmFilter "auto": percentile of the raw file's absolute mass errors (around its median) used as window
AUTO_PPM_FILTER_PERCENTILE = 90
# normalization "stable_drift": share of the precursors per retention time bin with the lowest spread over the runs
# that the drift is estimated from, width of the retention time bins (minutes) and minimum precursors per bin
STABLE_DRIFT_FRACTION = 0.75
STABLE_DRIFT_BIN_MINUTES = 0.1
STABLE_DRIFT_MIN_PER_BIN = 25
# driftBinMinutes "auto": stable_drift bin width as a share of the median peak width of the MS2 identifications (0.1 min
# was the best bin for 0.39 min wide peaks)
AUTO_DRIFT_BIN_PEAK_WIDTH = 0.25
# mbrRtWindow "auto": percentile of the retention time deviations of the MS2 identifications between the raw files
AUTO_MBR_RT_PERCENTILE = 90
# ppmTolerance "auto": mass error of a sample of identifications in the MS1 scan before their MS2 scan, searched
# within AUTO_PPM_SEARCH ppm (FlashLFQ's Orbitrap default) at the two most abundant isotopes, which have to agree within
# AUTO_PPM_ISOTOPE_AGREEMENT ppm; the tolerance is the |median| plus this percentile of the errors around the median
AUTO_PPM_SAMPLE = 5000
AUTO_PPM_SEARCH = 10.0
AUTO_PPM_ISOTOPE_AGREEMENT = 3.0
AUTO_PPM_PERCENTILE = 95
# FlashLFQ's defaults, used when there are too few identifications to measure the "auto" values
FLASHLFQ_PPM_TOLERANCE = 10.0
FLASHLFQ_MBR_RT_WINDOW = 1.5
MIN_AUTO_MEASUREMENTS = 100
PROTON = 1.007276467
C13_SPACING = 1.0033548
# natural abundances by nominal mass offset from the lightest isotope, and the compositions of common modifications
# (others are left out of the isotope envelope; their effect on which isotopes are the most abundant is small)
ISOTOPE_ABUNDANCES = {
    "C": [0.9893, 0.0107],
    "H": [0.999885, 0.000115],
    "N": [0.99636, 0.00364],
    "O": [0.99757, 0.00038, 0.00205],
    "S": [0.9499, 0.0075, 0.0425, 0.0, 0.0001],
    "P": [1.0],
}
UNIMOD_COMPOSITIONS = {
    "1": {"C": 2, "H": 2, "O": 1},  # acetyl
    "4": {"C": 2, "H": 3, "N": 1, "O": 1},  # carbamidomethyl
    "5": {"C": 1, "H": 1, "N": 1, "O": 1},  # carbamyl
    "7": {"H": -1, "N": -1, "O": 1},  # deamidation
    "21": {"H": 1, "O": 3, "P": 1},  # phospho
    "27": {"H": -2, "O": -1},  # Glu->pyro-Glu
    "28": {"H": -3, "N": -1},  # Gln->pyro-Glu
    "34": {"C": 1, "H": 2},  # methyl
    "35": {"O": 1},  # oxidation
    "36": {"C": 2, "H": 4},  # dimethyl
    "121": {"C": 4, "H": 6, "N": 2, "O": 2},  # GlyGly
}


def apply_quant(config: Config) -> None:
    """
    Quantify the rescored PSMs.

    With quantificationOptions.precursor "flashlfq", the MS1 signal of the rescored PSMs is measured with FlashLFQ,
    independent of the search engine (see :func:`quantify_ms1`). Otherwise ("search_engine"), the search engine's own
    quantification is updated with the rescored PSMs using picked-group-FDR (see :func:`quantify_search_engine`).

    :param config: config object containing all Oktoberfest parameters
    """
    if config.quant_precursor == "flashlfq":
        quantify_ms1(config)
    else:
        quantify_search_engine(config)


def quantify_search_engine(config: Config) -> None:
    """
    Call picked-group-FDR for means of quantification on search results with rescored Oktoberfest output.

    Integration of rescoring results only possible for MaxQuant at the moment. For msfragger and sage
    pciked_group_fdr is applied to the prior search results.

    :param config: config object containing all Oktoberfest parameters
    """
    if config.search_results.is_file():
        path = config.search_results.parent
    else:
        path = config.search_results

    fdr_dir = config.output / "results" / config.fdr_estimation_method

    # uses Percolator (or Mokapot) output as starting point
    okt_target_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.psms.txt")
    okt_decoy_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.decoy.psms.txt")

    # folder for picked group FDR output
    pg_fdr_dir = config.output / "picked_group_fdr/"
    pg_fdr_dir.mkdir(exist_ok=True, parents=True)

    fasta_file = str(config.library_input)

    digest_params_list = [
        DigestionParams(
            enzyme=config.enzyme,
            digestion=config.digestion,
            min_length=config.min_length,
            max_length=config.max_length,
            cleavages=config.missed_cleavages,
            special_aas=config.special_aas,
            fasta_contains_decoys=(True if config.db == "target" else False),
        )
    ]

    if config.search_results_type == "maxquant":
        evi_in = str(path / "evidence.txt")

        # saves updated MQ evidence file under same name, but different directory
        evi_out = str(pg_fdr_dir / "evidence.txt")

        logger.info("Update the MaxQuant evidence files with the rescored PSMs from Oktoberfest")
        pipeline.run_update_evidence(
            [evi_in], [okt_target_psms, okt_decoy_psms], [evi_out], "prosit", suppress_missing_peptide_warning=True
        )

        logger.info("Generate MQ proteinGroup-like file from updated evidence file")
        protein_groups_out = str(pg_fdr_dir / "rescore.proteinGroups.txt")

        lfq_min_peptide_ratios = 1  # 1 per default
        pipeline.run_picked_group_fdr(
            [evi_out],
            protein_groups_out,
            [fasta_file],
            digest_params_list,
            config.quantification,
            lfq_min_peptide_ratios,
            suppress_missing_peptide_warning=True,
        )

        logger.info("FDR-filtering the MQ proteinGroup-like file")
        fdr_cut = 0.01
        pipeline.run_filter_fdr_maxquant(
            [protein_groups_out],
            str(pg_fdr_dir / f"rescore.proteinGroups.fdr{int(fdr_cut * 100)}.txt"),
            fdr_cutoff=fdr_cut,
        )

    if config.search_results_type == "sage":
        cmd = [
            "python3", "-u", "-m", "picked_group_fdr",
            "--fasta", fasta_file,
            "--sage_results", str(path / "results.sage.tsv"),
            "--sage_lfq_tsv", str(path / "lfq.tsv"),
            "--protein_groups_out", str(pg_fdr_dir / "combined_protein.tsv"),
            "--output_format", "fragpipe",
            "--do_quant",
            "--lfq_min_peptide_ratios", "1",
            "--methods", "sage",
        ]  # fmt: skip
        logger.info(f"Starting picked-group-fdr with command {' '.join(cmd)}")
        subprocess.run(cmd, check=True)  # noqa: S603

    if config.search_results_type == "msfragger":
        # reusing IonQuant results by default
        cmd = [
            "python3", "-u", "-m", "picked_group_fdr",
            "--fasta", fasta_file,
            "--fragpipe_psm", *[str(p) for p in sorted(path.glob("*/psm.tsv"))],
            "--combined_ion", str(path / "combined_ion.tsv"),
            "--protein_groups_out", str(pg_fdr_dir / "combined_protein.tsv"),
            "--do_quant",
            "--lfq_min_peptide_ratios", "1",
            "--methods", "fragpipe",
        ]  # fmt: skip
        logger.info(f"Starting picked-group-fdr with command {' '.join(cmd)}")
        subprocess.run(cmd, check=True)  # noqa: S603


def quantify_ms1(config: Config) -> None:
    """
    Quantify the rescored PSMs on MS1, independent of the search engine.

    Writes to <output>/quantification:

    - precursors.tsv: one row per precursor (modified sequence + charge) with one normalized intensity column per raw
      file, in the ProteoBench "Custom" format (match-between-runs peaks only with quantificationOptions.mbrFdr)
    - peptidoforms.tsv: the precursors summed over their charge states, one row per modified sequence at the Percolator
      peptide q-value cutoff (quantificationOptions.peptideFdr); the precursors are only filtered by the PSM q-value
    - precursor_peaks.tsv: every FlashLFQ peak of an identified precursor, with its raw and normalized intensity,
      retention times, mass error, and whether it is used (or the reason why not)
    - normalization.tsv: the log2 correction curve of every raw file along the retention time ("drift",
      "stable_drift"), or the log2
      offset of every raw file to the first one ("mode", "median", "none")
    - protein_groups.tsv: directLFQ protein group intensities with the picked-group-FDR q-value
      (quantificationOptions.protein "directlfq")
    - settings.json: the FlashLFQ mass tolerance, match-between-runs retention time window and stable_drift bin width
      used, and the measurements behind the "auto" ones
    - flashlfq/, picked_group_fdr/, directlfq/: input and output of the tools

    :param config: config object containing all Oktoberfest parameters
    """
    quant_dir = config.output / "quantification"
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    targets = fdr_dir / "rescore.percolator.psms.txt"
    decoys = fdr_dir / "rescore.percolator.decoy.psms.txt"

    flashlfq_dir = quant_dir / "flashlfq"
    flashlfq_dir.mkdir(parents=True, exist_ok=True)
    ids = flashlfq_input(read_percolator_psms(targets), psm_metadata(config.output / "data"), config.quant_psm_fdr)
    ids.to_csv(flashlfq_dir / "ids.tsv", sep="\t", index=False)
    logger.info(f"Quantifying {ids['Full Sequence'].nunique()} precursors of {len(ids)} PSMs with FlashLFQ")
    link_spectra(config, ids["File Name"].unique(), flashlfq_dir / "spectra")
    settings = {
        "ppmTolerance": _setting(config.quant_ppm_tolerance, lambda: ppm_tolerance(ids, flashlfq_dir / "spectra",
                                                                                   config.num_threads)),
        "mbrRtWindow": _setting(config.quant_mbr_rt_window, lambda: mbr_rt_window(ids)),
    }  # fmt: skip
    run_flashlfq(config, flashlfq_dir, settings["ppmTolerance"]["value"], settings["mbrRtWindow"]["value"])

    peaks = pd.read_csv(flashlfq_dir / "QuantifiedPeaks.tsv", sep="\t", dtype={"File Name": str}, low_memory=False)
    runs = sorted(ids["File Name"].unique())
    peak_rows = peak_table(peaks, ids, config.quant_ppm_filter, config.quant_mbr_fdr)
    if config.quant_normalization == "stable_drift":
        settings["driftBinMinutes"] = _setting(config.quant_drift_bin_minutes, lambda: drift_bin_minutes(peak_rows))
    with open(quant_dir / "settings.json", "w") as f:
        json.dump(settings, f, indent=2)
    bin_minutes = settings.get("driftBinMinutes", {}).get("value", STABLE_DRIFT_BIN_MINUTES)
    normalization = normalize(peak_rows, runs, config.quant_normalization, bin_minutes)
    peak_rows.to_csv(quant_dir / "precursor_peaks.tsv", sep="\t", index=False)
    normalization.to_csv(quant_dir / "normalization.tsv", sep="\t", index=False)

    precursors = precursor_table(peak_rows, runs)
    precursors.to_csv(quant_dir / "precursors.tsv", sep="\t", index=False)
    peptides = peptides_at_fdr(fdr_dir / "rescore.percolator.peptides.txt", config.quant_peptide_fdr)
    peptidoform_table(precursors, peptides).to_csv(quant_dir / "peptidoforms.tsv", sep="\t", index=False)
    logger.info(f"Wrote {len(precursors)} quantified precursors to {quant_dir / 'precursors.tsv'}")

    picked_dir = quant_dir / "picked_group_fdr"
    protein_groups = run_picked_group_fdr(targets, decoys, config.library_input, picked_dir)
    protein_groups = protein_groups[(protein_groups["Q-value"] <= config.quant_protein_fdr)]
    protein_groups = protein_groups[protein_groups["Reverse"].isna()]
    logger.info(f"{len(protein_groups)} protein groups at {config.quant_protein_fdr} protein group FDR")

    if config.quant_protein == "directlfq":
        proteins = run_directlfq(precursors, protein_groups, quant_dir / "directlfq", config.num_threads)
        proteins.to_csv(quant_dir / "protein_groups.tsv", sep="\t", index=False)
        logger.info(f"Wrote {len(proteins)} protein group intensities to {quant_dir / 'protein_groups.tsv'}")


def read_percolator_psms(path: Path) -> pd.DataFrame:
    """Read a Percolator PSM table written by Oktoberfest (proteins joined by ";" in proteinIds)."""
    return pd.read_csv(path, sep="\t", dtype={"filename": str, "PSMId": str})


def psm_metadata(data_dir: Path) -> pd.DataFrame:
    """
    Read the retention time (minutes) and calculated mass of every PSM from the annotated spectra (hdf5).

    :param data_dir: <output>/data with one <raw file>.mzml.hdf5 per spectra file
    :return: one row per PSM with PSMId (as written to the Percolator tables), MODIFIED_SEQUENCE, PRECURSOR_CHARGE,
        RETENTION_TIME and CALCULATED_MASS
    """
    frames = []
    for path in sorted(data_dir.glob("*.mzml.hdf5")):
        obs = anndata.read_h5ad(path, backed="r").obs
        id_columns = SPEC_ID_COLUMNS + (["SCAN_EVENT_NUMBER"] if "SCAN_EVENT_NUMBER" in obs.columns else [])
        meta = obs[id_columns + ["RETENTION_TIME", "CALCULATED_MASS"]].copy()
        meta["PSMId"] = meta[id_columns].astype(str).agg("-".join, axis=1)
        frames.append(meta[["PSMId", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE", "RETENTION_TIME", "CALCULATED_MASS"]])
    if not frames:
        raise FileNotFoundError(f"No annotated spectra (*.mzml.hdf5) found in {data_dir}.")
    meta = pd.concat(frames, ignore_index=True).drop_duplicates("PSMId")
    meta["MODIFIED_SEQUENCE"] = meta["MODIFIED_SEQUENCE"].astype(str)
    return meta


def peptides_at_fdr(path: Path, fdr: float) -> set[str]:
    """Read the modified sequences at peptide q-value <= fdr from a Percolator peptide table ("_.PEPTIDE._")."""
    peptides = pd.read_csv(path, sep="\t", usecols=["peptide", "q-value"])
    return set(peptides.loc[peptides["q-value"] <= fdr, "peptide"].str.replace(r"^_\.|\._$", "", regex=True))


def target_proteins(protein_ids: str) -> str:
    """Remove the decoy proteins from a ";"-joined protein list."""
    return ";".join(p for p in str(protein_ids).split(";") if p and not p.startswith(DECOY_PREFIX))


def strip_modifications(modified_sequence: pd.Series) -> pd.Series:
    """Remove [UNIMOD:x] modifications (and the "-" of terminal ones) from modified sequences."""
    return modified_sequence.str.replace(r"\[[^\]]*\]-?", "", regex=True).str.replace("-", "", regex=False)


def flashlfq_input(psms: pd.DataFrame, meta: pd.DataFrame, psm_fdr: float) -> pd.DataFrame:
    """
    Build the FlashLFQ generic identification table from the target PSMs at q-value <= psm_fdr.

    FlashLFQ sums all charge states of one "Full Sequence", so the charge is appended to make it precursor specific.

    :param psms: Percolator target PSMs (PSMId, filename, q-value, proteinIds)
    :param meta: PSM metadata from :func:`psm_metadata`
    :param psm_fdr: PSM q-value cutoff
    :raises ValueError: if PSMs are missing in the annotated spectra
    :return: FlashLFQ identification table
    """
    psms = psms[psms["q-value"] <= psm_fdr].merge(meta, on="PSMId", how="left")
    missing = psms["RETENTION_TIME"].isna().sum()
    if missing:
        raise ValueError(f"{missing} rescored PSMs were not found in the annotated spectra (data/*.mzml.hdf5).")
    return pd.DataFrame(
        {
            "File Name": psms["filename"],
            "Base Sequence": strip_modifications(psms["MODIFIED_SEQUENCE"]),
            "Full Sequence": psms["MODIFIED_SEQUENCE"] + "/" + psms["PRECURSOR_CHARGE"].astype(int).astype(str),
            "Peptide Monoisotopic Mass": psms["CALCULATED_MASS"],
            "Scan Retention Time": psms["RETENTION_TIME"],
            "Precursor Charge": psms["PRECURSOR_CHARGE"].astype(int),
            "Protein Accession": psms["proteinIds"].map(target_proteins),
        }
    )


def link_spectra(config: Config, raw_files: list[str], spectra_dir: Path):
    """
    Link the mzML file of every raw file into one folder for FlashLFQ (which also writes its .ind index files there).

    mzML files converted by Oktoberfest are in <output>/spectra, given ones in the spectra input.

    :param config: config object containing all Oktoberfest parameters
    :param raw_files: raw file names (without extension) to quantify
    :param spectra_dir: folder to link the files into
    :raises FileNotFoundError: if no mzML file is found for a raw file
    """
    spectra_dir.mkdir(parents=True, exist_ok=True)
    search_dirs = [config.output / "spectra", config.spectra if config.spectra.is_dir() else config.spectra.parent]
    available = {p.stem: p for d in search_dirs if d.is_dir() for p in d.iterdir() if p.suffix.lower() == ".mzml"}
    for raw_file in raw_files:
        if raw_file not in available:
            raise FileNotFoundError(f"No mzML file with MS1 spectra found for {raw_file} in {search_dirs}.")
        link = spectra_dir / f"{raw_file}.mzML"
        if not link.exists():
            link.symlink_to(available[raw_file].resolve())


def _setting(configured: float | str, measure: Callable[[], dict]) -> dict:
    """Get a setting given in the config, or measured from the data ("auto"): {"value", "source", measurements}."""
    if configured == "auto":
        return {"source": "auto", **measure()}
    return {"value": float(configured), "source": "config"}


def ppm_tolerance(ids: pd.DataFrame, spectra_dir: Path, threads: int = 1) -> dict:
    """
    Measure FlashLFQ's MS1 mass tolerance from the MS1 spectra, before FlashLFQ runs.

    For a sample of AUTO_PPM_SAMPLE identifications, the MS1 scan before the MS2 scan is searched for the most intense
    peak within AUTO_PPM_SEARCH ppm of the most abundant isotope of the peptide whose second most abundant isotope is
    present at the same mass error (within AUTO_PPM_ISOTOPE_AGREEMENT ppm), as FlashLFQ needs an isotope envelope. The
    tolerance is the |median| error plus the AUTO_PPM_PERCENTILE percentile of the errors around the median, since
    FlashLFQ searches around the expected m/z.

    :param ids: FlashLFQ input from :func:`flashlfq_input`
    :param spectra_dir: folder with one <raw file>.mzML per raw file
    :param threads: number of raw files read in parallel
    :return: {"value": tolerance in ppm, "measured": identifications with a peak, "median_ppm", "percentile_ppm"};
        FlashLFQ's default FLASHLFQ_PPM_TOLERANCE if fewer than MIN_AUTO_MEASUREMENTS identifications have a peak
    """
    sample = ids.sample(n=min(AUTO_PPM_SAMPLE, len(ids)), random_state=42)
    groups = list(sample.groupby("File Name"))
    paths = [spectra_dir / f"{raw_file}.mzML" for raw_file, _ in groups]
    with ProcessPoolExecutor(max_workers=max(1, min(threads, len(groups)))) as executor:
        errors = [
            e for file_errors in executor.map(ms1_mass_errors, paths, [g for _, g in groups]) for e in file_errors
        ]
    return tolerance_from_errors(np.asarray(errors))


def tolerance_from_errors(errors: np.ndarray) -> dict:
    """Get the "auto" MS1 mass tolerance from measured mass errors (ppm), see :func:`ppm_tolerance`."""
    if len(errors) < MIN_AUTO_MEASUREMENTS:
        logger.warning(
            f"Only {len(errors)} MS1 mass errors measured, using FlashLFQ's default of {FLASHLFQ_PPM_TOLERANCE} ppm."
        )
        return {"value": FLASHLFQ_PPM_TOLERANCE, "measured": len(errors)}
    median = float(np.median(errors))
    spread = float(np.percentile(np.abs(errors - median), AUTO_PPM_PERCENTILE))
    value = round(abs(median) + spread, 2)
    logger.info(f"MS1 mass tolerance {value} ppm (median error {median:+.2f} ppm, {len(errors)} identifications)")
    return {"value": value, "measured": len(errors), "median_ppm": median, "percentile_ppm": spread}


def ms1_mass_errors(path: Path, ids: pd.DataFrame) -> list[float]:
    """
    Get the MS1 mass error (ppm) of identifications of one raw file, see :func:`ppm_tolerance`.

    :param path: mzML file
    :param ids: identifications of this raw file (Full Sequence, Peptide Monoisotopic Mass, Precursor Charge and
        Scan Retention Time in minutes)
    :return: mass error of every identification with an isotope envelope in the MS1 scan before its MS2 scan
    """
    with mzml.PreIndexedMzML(str(path), decode_binary=False) as reader:
        ms1 = [(float(s["scanList"]["scan"][0]["scan start time"]), s["id"]) for s in reader if s.get("ms level") == 1]
    ms1_rt = np.array([rt for rt, _ in ms1])
    errors = []
    columns = ["Scan Retention Time", "Full Sequence", "Peptide Monoisotopic Mass", "Precursor Charge"]
    with mzml.PreIndexedMzML(str(path)) as reader:
        for rt, full_sequence, mono_mass, charge in ids[columns].itertuples(index=False, name=None):
            k = np.searchsorted(ms1_rt, rt, side="right") - 1
            if k < 0:
                continue
            spectrum = reader.get_by_id(ms1[k][1])
            error = _envelope_error(
                spectrum["m/z array"], spectrum["intensity array"], full_sequence, mono_mass, int(charge)
            )
            if error is not None:
                errors.append(error)
    return errors


def _envelope_error(
    mz: np.ndarray, intensity: np.ndarray, full_sequence: str, mono_mass: float, charge: int
) -> Optional[float]:
    """Mass error of the most intense peak at the most abundant isotope with the second one at the same error."""
    first, second = np.argsort(-isotope_envelope(full_sequence))[:2]
    mono_mz = (mono_mass + charge * PROTON) / charge
    target = mono_mz + first * C13_SPACING / charge
    partner = mono_mz + second * C13_SPACING / charge
    ppm = (mz - target) / target * 1e6
    partner_ppm = (mz - partner) / partner * 1e6
    best, best_intensity = None, -1.0
    for i in np.flatnonzero(np.abs(ppm) <= AUTO_PPM_SEARCH):
        if intensity[i] > best_intensity and np.any(np.abs(partner_ppm - ppm[i]) <= AUTO_PPM_ISOTOPE_AGREEMENT):
            best, best_intensity = float(ppm[i]), intensity[i]
    return best


def isotope_envelope(full_sequence: str, peaks: int = 6) -> np.ndarray:
    """
    Get the relative abundances of the first isotopes of a peptide from its elemental composition.

    :param full_sequence: modified sequence in UNIMOD notation, optionally followed by "/<charge>"
    :param peaks: number of isotopes (0 = monoisotopic)
    :return: abundance of isotope 0..peaks-1
    """
    sequence = full_sequence.split("/")[0]
    composition = mass.Composition(sequence=re.sub(r"\[[^\]]*\]-?", "", sequence).replace("-", ""))
    for unimod in re.findall(r"\[UNIMOD:(\d+)\]", sequence):
        composition += mass.Composition(UNIMOD_COMPOSITIONS.get(unimod, {}))
    envelope = np.zeros(peaks)
    envelope[0] = 1.0
    for element, count in composition.items():
        for _ in range(max(count, 0)):
            envelope = np.convolve(envelope, ISOTOPE_ABUNDANCES.get(element, [1.0]))[:peaks]
    return envelope


def mbr_rt_window(ids: pd.DataFrame) -> dict:
    """
    Measure FlashLFQ's match-between-runs retention time window from the MS2 identifications.

    Every precursor identified in at least 3 raw files gets its median MS2 retention time per raw file. Its deviation
    from the precursor's median over the raw files, after removing every raw file's median deviation (a global shift),
    is how far a precursor elutes from where the other raw files predict it. The window is the AUTO_MBR_RT_PERCENTILE
    percentile of the absolute deviations.

    :param ids: FlashLFQ input from :func:`flashlfq_input`
    :return: {"value": window in minutes, "measured": deviations}; FlashLFQ's default FLASHLFQ_MBR_RT_WINDOW if there
        are fewer than MIN_AUTO_MEASUREMENTS deviations
    """
    rt = ids.groupby(["Full Sequence", "File Name"])["Scan Retention Time"].median().unstack()
    rt = rt[rt.notna().sum(axis=1) >= 3]
    deviation = rt.sub(rt.median(axis=1), axis=0)
    deviation = deviation.sub(deviation.median(), axis=1).abs().stack()
    if len(deviation) < MIN_AUTO_MEASUREMENTS:
        logger.warning(
            f"Only {len(deviation)} retention time deviations between raw files, using FlashLFQ's default "
            f"match-between-runs window of {FLASHLFQ_MBR_RT_WINDOW} min."
        )
        return {"value": FLASHLFQ_MBR_RT_WINDOW, "measured": len(deviation)}
    value = round(float(deviation.quantile(AUTO_MBR_RT_PERCENTILE / 100)), 3)
    logger.info(f"Match-between-runs retention time window {value} min ({len(deviation)} deviations)")
    return {"value": value, "measured": len(deviation)}


def drift_bin_minutes(peak_rows: pd.DataFrame) -> dict:
    """
    Get the stable_drift retention time bin from the peak width: AUTO_DRIFT_BIN_PEAK_WIDTH times the median width.

    :param peak_rows: peak table from :func:`peak_table`
    :return: {"value": bin width in minutes, "peak_width": median peak width of the used MS2 identifications}
    """
    used = peak_rows[peak_rows["Used"] & (peak_rows["Detection"] == "MSMS")]
    width = float((used["Peak RT end"] - used["Peak RT start"]).median())
    if not width > 0:
        logger.warning(f"No peak widths, using a stable_drift bin of {STABLE_DRIFT_BIN_MINUTES} min.")
        return {"value": STABLE_DRIFT_BIN_MINUTES}
    value = round(width * AUTO_DRIFT_BIN_PEAK_WIDTH, 4)
    logger.info(f"stable_drift bin {value} min (median peak width {width:.3f} min)")
    return {"value": value, "peak_width": width}


def run_flashlfq(
    config: Config,
    flashlfq_dir: Path,
    ppm: float = FLASHLFQ_PPM_TOLERANCE,
    mbr_rt_window: float = FLASHLFQ_MBR_RT_WINDOW,
):
    """
    Run FlashLFQ on <flashlfq_dir>/ids.tsv with the mzML files in <flashlfq_dir>/spectra.

    The FlashLFQ command line takes booleans as switches: --chg (quantify only the identified charge state) is on,
    normalization stays off, and match-between-runs (on by default, cannot be switched off) peaks are filtered later
    (see :func:`peak_table`). The isotope mass tolerance is the same as the peak one.

    :param config: config object containing all Oktoberfest parameters
    :param flashlfq_dir: folder with ids.tsv and spectra/, and for the output
    :param ppm: MS1 mass tolerance in ppm
    :param mbr_rt_window: maximum match-between-runs retention time window in minutes
    :raises RuntimeError: if FlashLFQ fails or writes no QuantifiedPeaks.tsv
    """
    exe = config.flashlfq_exe
    cmd = [config.dotnet_exe, str(exe)] if exe.suffix.lower() == ".dll" else [str(exe)]
    cmd += [
        "--idt", str(flashlfq_dir / "ids.tsv"),
        "--rep", str(flashlfq_dir / "spectra"),
        "--out", str(flashlfq_dir),
        "--chg",
        "--ppm", str(ppm),
        "--iso", str(ppm),
        "--mrt", str(mbr_rt_window),
        "--thr", str(config.num_threads),
    ]  # fmt: skip
    (flashlfq_dir / "QuantifiedPeaks.tsv").unlink(missing_ok=True)
    log_file = flashlfq_dir / "flashlfq.log"
    logger.info(f"Starting FlashLFQ with command {' '.join(cmd)}")
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)  # noqa: S603
    # FlashLFQ exits with 0 on most errors, so check its output
    if result.returncode != 0 or not (flashlfq_dir / "QuantifiedPeaks.tsv").is_file():
        raise RuntimeError(f"FlashLFQ failed (exit code {result.returncode}), see {log_file}.")


def peak_table(
    peaks: pd.DataFrame, ids: pd.DataFrame, ppm_filter: float | str | None = "auto", mbr_fdr: Optional[float] = None
) -> pd.DataFrame:
    """
    List every FlashLFQ peak of an identified precursor and decide which ones are used.

    Match-between-runs peaks are dropped, unless mbr_fdr is given: then those with a FlashLFQ PIP q-value <= mbr_fdr
    are kept, except FlashLFQ's decoy matches at a random retention time. A peak is used unless (in this order) it was
    not detected, it is shared by several precursors, its apex mass error is outside the window around the median of
    its raw file, or the precursor has a better peak in the same raw file (identifications far apart in retention
    time; a peak of an MS2 identification is preferred to a match-between-runs one, then the more intense one).
    The median mass error and the "auto" window of a raw file come from the peaks of its MS2 identifications.

    :param peaks: FlashLFQ QuantifiedPeaks.tsv
    :param ids: FlashLFQ input from :func:`flashlfq_input`
    :param ppm_filter: maximum distance of the apex mass error to the raw file's median in ppm; "auto" for
        the AUTO_PPM_FILTER_PERCENTILE percentile of the raw file's absolute mass errors around its median; None for no
        filter
    :param mbr_fdr: PIP q-value cutoff of the match-between-runs peaks, None to drop them
    :return: one row per peak and precursor with its "Detection" (MSMS or MBR), "Used" and the "Reason" if not
    """
    keep = peaks["Peak Detection Type"] == "MSMS"
    if mbr_fdr is not None:
        random_rt = peaks["Random RT"].astype(str).str.lower() == "true"
        q_value = pd.to_numeric(peaks["PIP Q-Value"], errors="coerce")
        keep |= (peaks["Peak Detection Type"] == "MBR") & ~random_rt & (q_value <= mbr_fdr)
    peaks = peaks[keep]
    rows = pd.DataFrame(
        {
            "Raw file": peaks["File Name"].astype(str),
            "Full Sequence": peaks["Full Sequence"].str.split("|", regex=False),
            "Intensity": pd.to_numeric(peaks["Peak intensity"], errors="coerce").fillna(0.0),
            "Detection": peaks["Peak Detection Type"],
            **{
                column: pd.to_numeric(peaks[source], errors="coerce")
                for column, source in [
                    ("Peak RT start", "Peak RT Start"),
                    ("Peak RT apex", "Peak RT Apex"),
                    ("Peak RT end", "Peak RT End"),
                    ("MS2 retention time", "MS2 Retention Time"),
                    ("Mass error (ppm)", "Peak Apex Mass Error (ppm)"),
                ]
            },
        }
    )
    rows["Shared"] = rows["Full Sequence"].str.len() > 1
    # median over the detected peaks before a shared peak is listed once per precursor
    run_median = (
        rows[(rows["Intensity"] > 0) & (rows["Detection"] == "MSMS")].groupby("Raw file")["Mass error (ppm)"].median()
    )
    rows["Mass error to run median (ppm)"] = rows["Mass error (ppm)"] - rows["Raw file"].map(run_median)
    rows["Mass error window (ppm)"] = rows["Raw file"].map(_mass_error_window(rows, ppm_filter))
    rows = rows.explode("Full Sequence", ignore_index=True)
    detected = rows["Intensity"] > 0

    reason = pd.Series("", index=rows.index)
    reason[~detected] = "not detected"
    reason[(reason == "") & rows["Shared"]] = "shared peak"
    if ppm_filter is not None:
        outside = ~(rows["Mass error to run median (ppm)"].abs() <= rows["Mass error window (ppm)"])
        reason[(reason == "") & outside] = "mass error"
    # "MBR" sorts before "MSMS", so the last peak is the most intense one of an MS2 identification if there is one
    candidates = rows[reason == ""].sort_values(["Detection", "Intensity"])
    best = candidates.groupby(["Full Sequence", "Raw file"]).tail(1).index
    reason[(reason == "") & ~rows.index.isin(best)] = "weaker peak"
    rows["Used"], rows["Reason"] = reason == "", reason

    precursor = rows["Full Sequence"].str.rsplit("/", n=1, expand=True)
    proteins = ids.groupby("Full Sequence")["Protein Accession"].first()
    rows.insert(1, "Modified sequence", precursor[0])
    rows.insert(2, "Charge", precursor[1].astype(int))
    rows.insert(3, "Proteins", rows["Full Sequence"].map(proteins))
    return rows.drop(columns=["Full Sequence", "Shared"]).sort_values(
        ["Modified sequence", "Charge", "Raw file"], ignore_index=True
    )


def _mass_error_window(rows: pd.DataFrame, ppm_filter: float | str | None) -> pd.Series:
    """Get the allowed distance of the apex mass error to the run median per raw file (NaN without a filter)."""
    runs = rows["Raw file"].unique()
    if ppm_filter is None:
        return pd.Series(np.nan, index=runs)
    if ppm_filter != "auto":
        return pd.Series(float(ppm_filter), index=runs)
    deviation = rows.loc[
        (rows["Intensity"] > 0) & (rows["Detection"] == "MSMS"), "Mass error to run median (ppm)"
    ].abs()
    window = deviation.groupby(rows["Raw file"]).quantile(AUTO_PPM_FILTER_PERCENTILE / 100).reindex(runs)
    for run, ppm in window.items():
        logger.info(f"Mass error window of {run}: +-{ppm:.2f} ppm around the median")
    return window


def normalize(
    peak_rows: pd.DataFrame, runs: list[str], normalization: str, bin_minutes: float = STABLE_DRIFT_BIN_MINUTES
) -> pd.DataFrame:
    """
    Add the "Normalized intensity" of every peak (in place), estimated from the MS2 identifications only.

    :param peak_rows: peak table from :func:`peak_table`
    :param runs: raw files
    :param normalization: "stable_drift", "drift", "mode", "median" or "none"
    :param bin_minutes: retention time bin of "stable_drift" in minutes
    :return: the normalization table: drift curves (:func:`drift_curves`) or offsets (:func:`normalization_table`)
    """
    identified = peak_rows[peak_rows["Detection"] == "MSMS"]
    if normalization in ("drift", "stable_drift"):
        if normalization == "stable_drift":
            curves = drift_curves(identified, runs, bin_minutes=bin_minutes, min_per_bin=STABLE_DRIFT_MIN_PER_BIN,
                                  stable_fraction=STABLE_DRIFT_FRACTION)  # fmt: skip
        else:
            curves = drift_curves(identified, runs)
        peak_rows["Normalized intensity"] = peak_rows["Intensity"] / np.power(2.0, drift_correction(peak_rows, curves))
        return curves
    offsets = normalization_offsets(identified, runs, normalization)
    peak_rows["Normalized intensity"] = peak_rows["Intensity"] / np.power(2.0, peak_rows["Raw file"].map(offsets))
    return normalization_table(identified, runs, offsets, normalization)


def normalization_offsets(peak_rows: pd.DataFrame, runs: list[str], normalization: str = "mode") -> pd.Series:
    """
    Get the log2 offset of every raw file to the first one, from the precursors used in both.

    "mode" takes the most frequent log2 ratio (histogram with 0.02 bins, smoothed over 9 bins), which follows the
    unchanged majority of the precursors; "median" takes the median log2 ratio, which is pulled by precursors that
    change between samples; "none" gives 0 for every raw file.

    :param peak_rows: peak table from :func:`peak_table`
    :param runs: raw files, the first one is the reference
    :param normalization: "mode", "median" or "none"
    :return: log2 offset per raw file (the intensities are divided by 2**offset)
    """
    used = peak_rows[peak_rows["Used"]]
    log = np.log2(
        used.pivot_table(index=["Modified sequence", "Charge"], columns="Raw file", values="Intensity", aggfunc="max")
    ).reindex(columns=runs)
    ratios = log.sub(log[runs[0]], axis=0)
    if normalization == "mode":
        offsets = ratios.apply(_mode)
    elif normalization == "median":
        offsets = ratios.median()
    else:
        offsets = pd.Series(0.0, index=runs)
    offsets[runs[0]] = 0.0
    return offsets.fillna(0.0)


def normalization_table(peak_rows: pd.DataFrame, runs: list[str], offsets: pd.Series, normalization: str):
    """
    Summarize the normalization of every raw file.

    One row per raw file with its log2 offset and the number of precursors it shares with the reference (first) raw
    file.

    :param peak_rows: peak table from :func:`peak_table`
    :param runs: raw files, the first one is the reference
    :param offsets: log2 offsets from :func:`normalization_offsets`
    :param normalization: normalization method
    :return: the normalization table
    """
    used = peak_rows[peak_rows["Used"]]
    found = used.pivot_table(index=["Modified sequence", "Charge"], columns="Raw file", values="Intensity",
                             aggfunc="max").reindex(columns=runs).notna()  # fmt: skip
    return pd.DataFrame(
        {
            "Raw file": runs,
            "Normalization": normalization,
            "Log2 offset": offsets.reindex(runs).values,
            "Precursors shared with reference": [int((found[runs[0]] & found[run]).sum()) for run in runs],
        }
    )


def drift_curves(
    peak_rows: pd.DataFrame,
    runs: list[str],
    min_runs: int = 3,
    bin_minutes: float = 0.5,
    min_per_bin: int = 100,
    stable_fraction: Optional[float] = None,
) -> pd.DataFrame:
    """
    Get the log2 intensity drift of every raw file along the retention time.

    Every precursor used in at least min_runs raw files (all of them if there are fewer) gets a reference level, its
    median log2 intensity over those raw files. Per raw file, the deviations from the reference are sorted by the
    apex retention time and grouped into bins of bin_minutes, each merged with its successors until it holds
    min_per_bin precursors. The curve is the deviation of every bin, smoothed over 3 bins. It holds the overall level
    of the raw file too, so no further normalization is needed and no raw file is the reference.

    Without stable_fraction ("drift"), all precursors are used and the deviation of a bin is its most frequent one
    (see :func:`_mode`), which follows the unchanged majority. With stable_fraction ("stable_drift"), only the most
    stable precursors are used (see :func:`stable_precursors`), so that precursors that change between the samples do
    not pull the curves, and the deviation of a bin is their median.

    :param peak_rows: peak table from :func:`peak_table`
    :param runs: raw files
    :param min_runs: minimum number of raw files a precursor is used in to take part
    :param bin_minutes: width of a retention time bin before merging, in minutes
    :param min_per_bin: minimum number of precursors in a bin
    :param stable_fraction: share of the precursors per retention time bin with the lowest spread over the raw files
        that the curves are estimated from, None to use all precursors
    :return: one row per raw file and bin: Raw file, Normalization, Retention time (median of the bin), Log2 offset
        (the curve) and Precursors (in the bin)
    """
    used = peak_rows[peak_rows["Used"]]
    index = ["Modified sequence", "Charge"]
    log = np.log2(used.pivot_table(index=index, columns="Raw file", values="Intensity", aggfunc="max")).reindex(
        columns=runs
    )
    rt = used.pivot_table(index=index, columns="Raw file", values="Peak RT apex", aggfunc="median").reindex(
        index=log.index, columns=runs
    )
    usable = log.notna().sum(axis=1) >= min(min_runs, len(runs))
    log, rt = log[usable], rt[usable]
    if stable_fraction is None:
        return _drift_curves(log, rt, runs, bin_minutes, min_per_bin, "drift", _mode)
    stable = stable_precursors(log, rt, stable_fraction, bin_minutes, min_per_bin)
    return _drift_curves(log[stable], rt[stable], runs, bin_minutes, min_per_bin, "stable_drift", _median)


def _drift_curves(
    log: pd.DataFrame,
    rt: pd.DataFrame,
    runs: list[str],
    bin_minutes: float,
    min_per_bin: int,
    name: str,
    center: Callable[[pd.Series], float],
) -> pd.DataFrame:
    """Drift curves of the precursors (rows) of log (log2 intensities) and rt (apex retention times), see drift_curves."""
    deviations = log.sub(log.median(axis=1), axis=0)
    curves = []
    for run in runs:
        points = pd.DataFrame({"rt": rt[run], "deviation": deviations[run]}).dropna().sort_values("rt")
        if points.empty:
            continue
        points["bin"] = _merge_bins(np.floor(points["rt"].to_numpy() / bin_minutes), min_per_bin)
        bins = points.groupby("bin").agg(rt=("rt", "median"), precursors=("rt", "size"))
        bins["offset"] = points.groupby("bin")["deviation"].apply(center)
        curves.append(
            pd.DataFrame(
                {
                    "Raw file": run,
                    "Normalization": name,
                    "Retention time": bins["rt"].to_numpy(),
                    "Log2 offset": bins["offset"].rolling(3, center=True, min_periods=1).mean().to_numpy(),
                    "Precursors": bins["precursors"].to_numpy(),
                }
            )
        )
    columns = ["Raw file", "Normalization", "Retention time", "Log2 offset", "Precursors"]
    return pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(columns=columns)


def stable_precursors(
    log: pd.DataFrame,
    rt: pd.DataFrame,
    fraction: float,
    bin_minutes: float = 0.1,
    min_per_bin: int = 25,
) -> pd.Series:
    """
    Select the precursors whose intensity is most stable over the raw files.

    Every precursor gets the standard deviation of its log2 intensities over the raw files. Precursors of proteins that
    change between the samples have a larger spread, without the samples having to be known. The precursors are
    grouped by their median apex retention time (bins of bin_minutes merged until they hold min_per_bin / fraction
    precursors), and the fraction of every bin with the lowest standard deviation is kept, so that every part of the
    gradient keeps precursors.

    :param log: log2 intensities, one row per precursor and one column per raw file
    :param rt: apex retention times like log
    :param fraction: share of the precursors kept per retention time bin
    :param bin_minutes: width of a retention time bin before merging, in minutes
    :param min_per_bin: minimum number of kept precursors in a bin
    :return: boolean mask of the stable precursors (index of log)
    """
    spread = log.std(axis=1)
    center = rt.median(axis=1)
    has_spread = spread.notna() & center.notna()
    bins = pd.Series(np.nan, index=log.index)
    bins[has_spread] = _merge_bins(
        np.floor(center[has_spread].to_numpy() / bin_minutes), int(np.ceil(min_per_bin / fraction))
    )
    rank = spread[has_spread].groupby(bins[has_spread]).rank(pct=True, method="first")
    return (rank <= fraction).reindex(log.index, fill_value=False)


def _merge_bins(fine_bins: np.ndarray, min_per_bin: int) -> np.ndarray:
    """Merge sorted fine bin numbers into bins of at least min_per_bin entries; a short last bin joins the one before."""
    values, counts = np.unique(fine_bins, return_counts=True)
    merged, current, count = np.zeros(len(values), dtype=int), 0, 0
    for i, n in enumerate(counts):
        merged[i] = current
        count += n
        if count >= min_per_bin:
            current, count = current + 1, 0
    if count and current:
        merged[merged == current] = current - 1
    return merged[np.searchsorted(values, fine_bins)]


def drift_correction(peak_rows: pd.DataFrame, curves: pd.DataFrame) -> pd.Series:
    """
    Get the log2 correction of every peak: the drift curve of its raw file at its apex retention time.

    :param peak_rows: peak table from :func:`peak_table`
    :param curves: drift curves from :func:`drift_curves`
    :return: log2 correction per peak row (the intensities are divided by 2**correction), 0 without a curve or apex
    """
    correction = pd.Series(0.0, index=peak_rows.index)
    for run, curve in curves.groupby("Raw file"):
        rows = peak_rows["Raw file"] == run
        correction[rows] = np.interp(peak_rows.loc[rows, "Peak RT apex"], curve["Retention time"], curve["Log2 offset"])
    return correction.fillna(0.0)


def _median(log_ratios: pd.Series) -> float:
    """Median of a log2 ratio distribution, 0 without values."""
    log_ratios = log_ratios.dropna()
    return float(log_ratios.median()) if not log_ratios.empty else 0.0


def _mode(log_ratios: pd.Series) -> float:
    """
    Most frequent value of a log2 ratio distribution.

    The highest bin of a histogram (0.02 bins in [-3, 3], smoothed over 9 bins) is refined to the median of the ratios
    within 0.1 of it.
    """
    log_ratios = log_ratios.dropna()
    if log_ratios.empty:
        return 0.0
    hist, edges = np.histogram(log_ratios, bins=np.arange(-3, 3.0001, 0.02))
    smoothed = np.convolve(hist, np.ones(9) / 9, mode="same")
    # middle of the highest bins, so that a flat top is not resolved to its left edge
    i = int(round(np.flatnonzero(np.isclose(smoothed, smoothed.max())).mean()))
    center = (edges[i] + edges[i + 1]) / 2
    return float(log_ratios[(log_ratios - center).abs() <= 0.1].median())


def precursor_table(peak_rows: pd.DataFrame, runs: list[str]) -> pd.DataFrame:
    """
    Collapse the used peaks to one normalized intensity per precursor and raw file.

    :param peak_rows: peak table from :func:`peak_table` with "Normalized intensity"
    :param runs: raw files (intensity columns)
    :return: one row per precursor: Sequence, Modified sequence, Proteins, Charge and one intensity column per raw file
        (the ProteoBench "Custom" format)
    """
    used = peak_rows[peak_rows["Used"]]
    intensities = used.pivot_table(
        index=["Modified sequence", "Charge"], columns="Raw file", values="Normalized intensity", aggfunc="max"
    ).reindex(columns=runs)
    proteins = used.groupby(["Modified sequence", "Charge"])["Proteins"].first()
    table = intensities.index.to_frame(index=False)
    table.insert(0, "Sequence", strip_modifications(table["Modified sequence"]))
    table.insert(2, "Proteins", proteins.reindex(intensities.index).values)
    for run in runs:
        table[run] = intensities[run].values
    return table


def peptidoform_table(precursors: pd.DataFrame, peptides: set[str] | None = None) -> pd.DataFrame:
    """
    Sum the precursor intensities of every modified sequence over its charge states.

    :param precursors: precursor table from :func:`precursor_table`
    :param peptides: if given, only these modified sequences (e.g. those at the peptide q-value cutoff)
    :return: one row per modified sequence: Sequence, Modified sequence, Proteins and one intensity column per raw file
        (the ProteoBench "Custom" format of the peptidoform modules)
    """
    runs = list(precursors.columns[4:])
    if peptides is not None:
        precursors = precursors[precursors["Modified sequence"].isin(peptides)]
    grouped = precursors.groupby("Modified sequence", sort=True)
    table = grouped[["Sequence", "Proteins"]].first().reset_index()[["Sequence", "Modified sequence", "Proteins"]]
    sums = grouped[runs].sum(min_count=1)
    for run in runs:
        table[run] = sums[run].values
    return table


def run_picked_group_fdr(targets: Path, decoys: Path, fasta: Path, out_dir: Path) -> pd.DataFrame:
    """
    Infer protein groups with picked-group-FDR from the Percolator target and decoy PSMs.

    picked-group-FDR expects one protein per trailing column, Oktoberfest writes them joined by ";" in proteinIds.

    :param targets: Percolator target PSMs
    :param decoys: Percolator decoy PSMs
    :param fasta: fasta file of the search
    :param out_dir: output folder
    :return: the protein groups (proteinGroups.txt) with Q-value and Reverse columns
    """
    import picked_group_fdr.picked_group_fdr as picked

    out_dir.mkdir(parents=True, exist_ok=True)
    evidence = [split_protein_columns(path, out_dir / path.name) for path in (targets, decoys)]
    protein_groups = out_dir / "proteinGroups.txt"
    picked.main(
        ["--perc_evidence", *map(str, evidence), "--methods", "picked_protein_group_no_remap",
         "--fasta", str(fasta), "--protein_groups_out", str(protein_groups)]
    )  # fmt: skip
    return pd.read_csv(protein_groups, sep="\t")


def split_protein_columns(src: Path, dst: Path) -> Path:
    """Write a Percolator PSM table with the ";"-joined proteinIds split into one protein per trailing column."""
    with open(src) as fin, open(dst, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        protein_column = header.index("proteinIds")
        fout.write("\t".join(header) + "\n")
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            fields = fields[:protein_column] + fields[protein_column].split(";")
            fout.write("\t".join(fields) + "\n")
    return dst


def precursors_by_protein_group(precursors: pd.DataFrame, protein_groups: pd.DataFrame) -> pd.DataFrame:
    """
    Assign every precursor to its protein group; precursors of several or no accepted protein groups are dropped.

    :param precursors: precursor table from :func:`precursor_table`
    :param protein_groups: accepted protein groups ("Protein IDs" joined by ";")
    :return: directLFQ input: protein, ion and one intensity column per raw file (missing values as 0)
    """
    protein_to_group = {p: group for group in protein_groups["Protein IDs"] for p in group.split(";")}

    def group_of(proteins: str):
        groups = {protein_to_group.get(p) for p in str(proteins).split(";")} - {None}
        return groups.pop() if len(groups) == 1 else None

    runs = list(precursors.columns[4:])
    table = precursors.assign(
        protein=precursors["Proteins"].map(group_of),
        ion=precursors["Modified sequence"] + "/" + precursors["Charge"].astype(str),
    ).dropna(subset=["protein"])
    logger.info(f"{len(table)} of {len(precursors)} precursors belong to exactly one accepted protein group")
    return table[["protein", "ion"] + runs].fillna(0)


def run_directlfq(precursors: pd.DataFrame, protein_groups: pd.DataFrame, out_dir: Path, threads: int) -> pd.DataFrame:
    """
    Calculate protein group intensities with directLFQ.

    :param precursors: precursor table from :func:`precursor_table`
    :param protein_groups: accepted protein groups
    :param out_dir: folder for the directLFQ input and output
    :param threads: number of threads
    :raises ImportError: if directLFQ is not installed
    :return: one row per protein group: Protein group, Q-value, Precursors (number used) and one intensity column per
        raw file
    """
    try:
        from directlfq import lfq_manager
    except ImportError as e:
        raise ImportError(
            "Protein quantification with directLFQ needs the directlfq package: pip install oktoberfest[quant], or "
            'set quantificationOptions.protein to "none".'
        ) from e

    out_dir.mkdir(parents=True, exist_ok=True)
    # the .aq_reformat.tsv suffix tells directLFQ the table is in its own format already
    aq_file = out_dir / "precursors.aq_reformat.tsv"
    aq = precursors_by_protein_group(precursors, protein_groups)
    aq.to_csv(aq_file, sep="\t", index=False)
    lfq_manager.run_lfq(str(aq_file), num_cores=threads)
    proteins = pd.read_csv(f"{aq_file}.protein_intensities.tsv", sep="\t")
    runs = [column for column in proteins.columns if column != "protein"]
    proteins.insert(1, "Q-value", proteins["protein"].map(protein_groups.set_index("Protein IDs")["Q-value"]))
    proteins.insert(2, "Precursors", proteins["protein"].map(aq["protein"].value_counts()).astype(int))
    proteins[runs] = proteins[runs].replace(0, np.nan)
    return proteins.rename(columns={"protein": "Protein group"})
