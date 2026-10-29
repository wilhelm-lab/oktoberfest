import logging
import subprocess
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from picked_group_fdr.digestion_params import DigestionParams
from picked_group_fdr.pipeline import pipeline

from . import Config

logger = logging.getLogger(__name__)

DECOY_PREFIX = "rev_"
SPEC_ID_COLUMNS = ["RAW_FILE", "SCAN_NUMBER", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"]


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
      file, in the ProteoBench "Custom" format (no match-between-runs)
    - peptidoforms.tsv: the precursors summed over their charge states, one row per modified sequence at the Percolator
      peptide q-value cutoff (quantificationOptions.peptideFdr); the precursors are only filtered by the PSM q-value
    - precursor_peaks.tsv: every FlashLFQ peak of an identified precursor, with its raw and normalized intensity,
      retention times, mass error, and whether it is used (or the reason why not)
    - normalization.tsv: the log2 correction curve of every raw file along the retention time ("drift"), or the log2
      offset of every raw file to the first one ("mode", "median", "none")
    - protein_groups.tsv: directLFQ protein group intensities with the picked-group-FDR q-value
      (quantificationOptions.protein "directlfq")
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
    run_flashlfq(config, flashlfq_dir)

    peaks = pd.read_csv(flashlfq_dir / "QuantifiedPeaks.tsv", sep="\t", dtype={"File Name": str}, low_memory=False)
    runs = sorted(ids["File Name"].unique())
    peak_rows = peak_table(peaks, ids, config.quant_ppm_filter)
    if config.quant_normalization == "drift":
        curves = drift_curves(peak_rows, runs)
        peak_rows["Normalized intensity"] = peak_rows["Intensity"] / np.power(2.0, drift_correction(peak_rows, curves))
        normalization = curves
    else:
        offsets = normalization_offsets(peak_rows, runs, config.quant_normalization)
        peak_rows["Normalized intensity"] = peak_rows["Intensity"] / np.power(2.0, peak_rows["Raw file"].map(offsets))
        normalization = normalization_table(peak_rows, runs, offsets, config.quant_normalization)
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


def run_flashlfq(config: Config, flashlfq_dir: Path):
    """
    Run FlashLFQ on <flashlfq_dir>/ids.tsv with the mzML files in <flashlfq_dir>/spectra.

    The FlashLFQ command line takes booleans as switches: --chg (quantify only the identified charge state) is on,
    normalization stays off, and match-between-runs (on by default, cannot be switched off) peaks are dropped later.

    :param config: config object containing all Oktoberfest parameters
    :param flashlfq_dir: folder with ids.tsv and spectra/, and for the output
    :raises RuntimeError: if FlashLFQ fails or writes no QuantifiedPeaks.tsv
    """
    exe = config.flashlfq_exe
    cmd = [config.dotnet_exe, str(exe)] if exe.suffix.lower() == ".dll" else [str(exe)]
    cmd += [
        "--idt", str(flashlfq_dir / "ids.tsv"),
        "--rep", str(flashlfq_dir / "spectra"),
        "--out", str(flashlfq_dir),
        "--chg",
        "--ppm", str(config.quant_ppm_tolerance),
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


def peak_table(peaks: pd.DataFrame, ids: pd.DataFrame, ppm_filter: float | None = 5) -> pd.DataFrame:
    """
    List every FlashLFQ peak of an identified precursor and decide which ones are used.

    Match-between-runs peaks are dropped. A peak is used unless (in this order) it was not detected, it is shared by
    several precursors, its apex mass error is more than ppm_filter from the median of its raw file, or the precursor
    has a more intense peak in the same raw file (identifications far apart in retention time).

    :param peaks: FlashLFQ QuantifiedPeaks.tsv
    :param ids: FlashLFQ input from :func:`flashlfq_input`
    :param ppm_filter: maximum distance of the apex mass error to the raw file's median in ppm, None for no filter
    :return: one row per peak and precursor with "Used" and the "Reason" if not
    """
    peaks = peaks[peaks["Peak Detection Type"] == "MSMS"]
    rows = pd.DataFrame(
        {
            "Raw file": peaks["File Name"].astype(str),
            "Full Sequence": peaks["Full Sequence"].str.split("|", regex=False),
            "Intensity": pd.to_numeric(peaks["Peak intensity"], errors="coerce").fillna(0.0),
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
    run_median = rows[rows["Intensity"] > 0].groupby("Raw file")["Mass error (ppm)"].median()
    rows["Mass error to run median (ppm)"] = rows["Mass error (ppm)"] - rows["Raw file"].map(run_median)
    rows = rows.explode("Full Sequence", ignore_index=True)
    detected = rows["Intensity"] > 0

    reason = pd.Series("", index=rows.index)
    reason[~detected] = "not detected"
    reason[(reason == "") & rows["Shared"]] = "shared peak"
    if ppm_filter is not None:
        reason[(reason == "") & ~(rows["Mass error to run median (ppm)"].abs() <= ppm_filter)] = "mass error"
    candidates = rows[reason == ""]
    strongest = candidates.groupby(["Full Sequence", "Raw file"])["Intensity"].idxmax()
    reason[(reason == "") & ~rows.index.isin(strongest)] = "weaker peak"
    rows["Used"], rows["Reason"] = reason == "", reason

    precursor = rows["Full Sequence"].str.rsplit("/", n=1, expand=True)
    proteins = ids.groupby("Full Sequence")["Protein Accession"].first()
    rows.insert(1, "Modified sequence", precursor[0])
    rows.insert(2, "Charge", precursor[1].astype(int))
    rows.insert(3, "Proteins", rows["Full Sequence"].map(proteins))
    return rows.drop(columns=["Full Sequence", "Shared"]).sort_values(
        ["Modified sequence", "Charge", "Raw file"], ignore_index=True
    )


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
    peak_rows: pd.DataFrame, runs: list[str], min_runs: int = 3, bin_minutes: float = 0.5, min_per_bin: int = 100
) -> pd.DataFrame:
    """
    Get the log2 intensity drift of every raw file along the retention time.

    Every precursor used in at least min_runs raw files (all of them if there are fewer) gets a reference level, its
    median log2 intensity over those raw files. Per raw file, the deviations from the reference are sorted by the
    apex retention time and grouped into bins of bin_minutes, each merged with its successors until it holds
    min_per_bin precursors. The curve is the most frequent deviation of every bin (see :func:`_mode`), smoothed over 3
    bins. It holds the overall level of the raw file too, so no further normalization is needed and no raw file is the
    reference.

    :param peak_rows: peak table from :func:`peak_table`
    :param runs: raw files
    :param min_runs: minimum number of raw files a precursor is used in to take part
    :param bin_minutes: width of a retention time bin before merging, in minutes
    :param min_per_bin: minimum number of precursors in a bin
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
    deviations = log[usable].sub(log[usable].median(axis=1), axis=0)

    curves = []
    for run in runs:
        points = pd.DataFrame({"rt": rt.loc[usable, run], "deviation": deviations[run]}).dropna().sort_values("rt")
        if points.empty:
            continue
        points["bin"] = _merge_bins(np.floor(points["rt"].to_numpy() / bin_minutes), min_per_bin)
        bins = points.groupby("bin").agg(rt=("rt", "median"), precursors=("rt", "size"))
        bins["offset"] = points.groupby("bin")["deviation"].apply(_mode)
        curves.append(
            pd.DataFrame(
                {
                    "Raw file": run,
                    "Normalization": "drift",
                    "Retention time": bins["rt"].to_numpy(),
                    "Log2 offset": bins["offset"].rolling(3, center=True, min_periods=1).mean().to_numpy(),
                    "Precursors": bins["precursors"].to_numpy(),
                }
            )
        )
    columns = ["Raw file", "Normalization", "Retention time", "Log2 offset", "Precursors"]
    return pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(columns=columns)


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
