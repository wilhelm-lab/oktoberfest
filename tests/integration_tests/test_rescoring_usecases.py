"""
Integration tests for the Oktoberfest rescoring use cases.

Every test runs only an Oktoberfest Rescoring job on prepared inputs: MS2-only mzML files (ThermoRawFileParser) and
MSFragger 4.2 search results (pepXML), which are downloaded from figshare once and cached in
tests/data/<usecase>/inputs. The checks are deliberately robust (output files, feature columns, more PSMs at 1% FDR
after rescoring than before), since Koina predictions and Percolator vary slightly between runs.

The download URL of a use case can be overridden with the environment variable OKTOBERFEST_TEST_URL_<USECASE>
(e.g. OKTOBERFEST_TEST_URL_TMT), which accepts http(s):// and file:// URLs as well as local paths to the zip file.
"""

import json
import os
import re
import shutil
import unittest
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from oktoberfest.__main__ import main

from .. import CONFIGS_DIR, DATA_DIR

# figshare item 10.6084/m9.figshare.34876347 (one zip per use case)
FIGSHARE_URL_TMT = "https://ndownloader.figshare.com/files/70952373"
FIGSHARE_URL_CIT = "https://ndownloader.figshare.com/files/70952382"
FIGSHARE_URL_MULTIFRAG = "https://ndownloader.figshare.com/files/70952376"
FIGSHARE_URL_PHOSPHO = "https://ndownloader.figshare.com/files/70952370"
FIGSHARE_URL_HLA = "https://ndownloader.figshare.com/files/70952379"
FIGSHARE_URL_QUANT = "https://ndownloader.figshare.com/files/71027736"

FIGSHARE_URLS = {
    "tmt": FIGSHARE_URL_TMT,
    "cit": FIGSHARE_URL_CIT,
    "multifrag": FIGSHARE_URL_MULTIFRAG,
    "phospho": FIGSHARE_URL_PHOSPHO,
    "hla": FIGSHARE_URL_HLA,
    "quant": FIGSHARE_URL_QUANT,
}

FDR = 0.01


def _download_inputs(usecase: str) -> Path:
    """
    Download and extract the inputs of a use case into tests/data/<usecase>/inputs, unless they are already there.

    :param usecase: name of the use case
    :raises unittest.SkipTest: if no download URL is configured yet
    :return: the inputs directory
    """
    inputs_dir = DATA_DIR / usecase / "inputs"
    done_marker = inputs_dir / ".download_complete"
    if done_marker.is_file():
        return inputs_dir
    url = os.environ.get(f"OKTOBERFEST_TEST_URL_{usecase.upper()}", FIGSHARE_URLS[usecase])
    if url.startswith("TODO"):
        raise unittest.SkipTest(f"No download URL for the {usecase} inputs configured yet ({url}).")
    inputs_dir.mkdir(parents=True, exist_ok=True)
    zip_path = inputs_dir / f"{usecase}.zip"
    scheme = urllib.parse.urlparse(url).scheme
    if scheme in ("http", "https"):
        urllib.request.urlretrieve(url, zip_path)  # noqa: S310 (only http and https URLs are opened)
    else:
        # local zip file, given as a file:// URL or a path
        local_path = urllib.request.url2pathname(urllib.parse.urlparse(url).path) if scheme == "file" else url
        shutil.copy(local_path, zip_path)
    with zipfile.ZipFile(zip_path) as zip_file:
        zip_file.extractall(inputs_dir)
    zip_path.unlink()
    done_marker.touch()
    return inputs_dir


def _run_rescoring(config_name: str):
    """Run the rescoring job of a test config through the command line entry point."""
    config_path = CONFIGS_DIR / config_name
    with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
        main()


def _accepted(psms_file: Path) -> pd.DataFrame:
    """Read a percolator target PSM or peptide table and keep the entries at 1% FDR."""
    df = pd.read_csv(psms_file, sep="\t", dtype={"filename": str})
    return df[df["q-value"] < FDR]


def _n_accepted(*psms_files: Path) -> int:
    """Count the PSMs (or peptides) at 1% FDR over one or more percolator target tables."""
    return sum(len(_accepted(psms_file)) for psms_file in psms_files)


def _tab_columns(tab_file: Path) -> list[str]:
    """Read the header of a percolator tab file."""
    with open(tab_file) as f:
        return f.readline().rstrip("\n").split("\t")


def _species_ratios(table: pd.DataFrame, condition_a: list[str], condition_b: list[str], expected: dict) -> dict:
    """
    Get the log2 A/B ratios per species and the absolute error to the expected ratio, as ProteoBench computes them.

    Entries quantified in at least 3 runs with a value in both conditions; the ratio is the log2 of the mean intensity
    in A over the mean in B. Entries of more than one species are left out.

    :param table: one row per precursor or protein, a "Proteins" or "Protein group" column and one column per run
    :param condition_a: runs of condition A
    :param condition_b: runs of condition B
    :param expected: expected log2 A/B ratio per species (HUMAN, YEAST, ECOLI)
    :return: number of entries, median log2 ratio per species and median absolute error (global and per species)
    """
    proteins = table["Proteins"] if "Proteins" in table else table["Protein group"]
    species = proteins.map(lambda p: {name.rsplit("_", 1)[-1] for name in str(p).split(";")})
    species = species.map(lambda s: next(iter(s)) if len(s) == 1 else None)
    runs = condition_a + condition_b
    used = (
        (table[runs].notna().sum(axis=1) >= 3)
        & table[condition_a].notna().any(axis=1)
        & table[condition_b].notna().any(axis=1)
        & species.isin(list(expected))
    )
    log2 = np.log2(table.loc[used, condition_a].mean(axis=1) / table.loc[used, condition_b].mean(axis=1))
    epsilon = (log2 - species[used].map(expected)).abs()
    result = {"entries": int(used.sum()), "median_abs_epsilon": float(epsilon.median())}
    for name in expected:
        result[f"log2_{name}"] = float(log2[species[used] == name].median())
        result[f"median_abs_epsilon_{name}"] = float(epsilon[species[used] == name].median())
    return result


def _phospho_sites(peptide: str) -> tuple[str, set[int]]:
    """Get the plain sequence and the 1-based phospho sites of a percolator peptide like _.AS[UNIMOD:21]PEK._."""
    sequence, sites = "", set()
    peptide = re.sub(r"^\[UNIMOD:\d+\]-", "", peptide.strip("_."))
    for residue, mod in re.findall(r"([A-Z])((?:\[UNIMOD:\d+\])*)", peptide):
        sequence += residue
        if "UNIMOD:21" in mod:
            sites.add(len(sequence))
    return sequence, sites


class TestRescoringUseCases(unittest.TestCase):
    """Rescoring use cases on MSFragger results: TMT, citrullination, MultiFrag, phospho localization, HLA PTM FDR."""

    def _prepare(self, usecase: str, out_dirs: list[str]):
        """Download the inputs and register the clean-up of the outputs."""
        _download_inputs(usecase)
        for out_dir in out_dirs:
            out_path = DATA_DIR / usecase / out_dir
            shutil.rmtree(out_path, ignore_errors=True)
            self.addCleanup(shutil.rmtree, out_path, ignore_errors=True)

    @staticmethod
    def _expected(usecase: str, name: str = "expected.json") -> dict:
        """Load the expected outputs of a use case."""
        with open(DATA_DIR / usecase / "expected_outputs" / name) as f:
            return json.load(f)

    def _check_files(self, out_dir: Path, files: list[str]):
        for file in files:
            self.assertTrue((out_dir / file).is_file(), f"Missing output file {out_dir / file}")

    def _check_columns(self, tab_file: Path, expected_columns: list[str]):
        self.assertListEqual(_tab_columns(tab_file), expected_columns, f"Unexpected columns in {tab_file}")

    def _check_range(self, value: float, value_range: list[float], label: str):
        low, high = value_range
        self.assertTrue(low <= value <= high, f"{label}: {value} outside [{low}, {high}]")

    def _check_gain(self, original: int, rescored: int, label: str):
        self.assertGreater(original, 0, f"No {label} at 1% FDR before rescoring")
        self.assertGreater(rescored, original, f"Rescoring lost {label} at 1% FDR: {original} -> {rescored}")

    def _check_standard_run(self, out_dir: Path, expected: dict):
        """Check the outputs of a run without PTM-specific FDR control or localization."""
        fdr_dir = out_dir / "results" / "percolator"
        self._check_files(out_dir, expected["output_files"])
        self._check_columns(fdr_dir / "original.tab", expected["original_columns"])
        self._check_columns(fdr_dir / "rescore.tab", expected["rescore_columns"])
        for level in ("psms", "peptides"):
            self._check_gain(
                _n_accepted(fdr_dir / f"original.percolator.{level}.txt"),
                _n_accepted(fdr_dir / f"rescore.percolator.{level}.txt"),
                level,
            )
        return fdr_dir

    def test_rescoring_tmt(self):
        """Rescoring of TMT10 labelled peptides (ion trap CID) with the Prosit TMT models (PXD030339)."""
        self._prepare("tmt", ["out"])
        expected = self._expected("tmt")
        _run_rescoring("rescoring_tmt.json")
        fdr_dir = self._check_standard_run(DATA_DIR / "tmt" / "out", expected)
        # all accepted PSMs carry the TMT label on the N-terminus
        peptides = _accepted(fdr_dir / "rescore.percolator.psms.txt")["peptide"]
        self.assertTrue(peptides.str.contains(r"\[UNIMOD:737\]-", regex=True).all())

    def test_rescoring_citrullination(self):
        """
        Citrullination localization (UNIMOD:7 on N, Q and R) with the Prosit Cit models (PXD056560).

        Synthetic citrullinated peptides spiked into a MES-SA digest. The citrullination workflow of localize_ptms
        keeps only spectra with a UNIMOD:7 candidate and rescores them in localize_mod. original.tab contains all
        localization candidates, so the original percolator results are no valid baseline; the number of
        citrullinated PSMs is compared with the reference run instead (the PRIDE deposit of the Prosit-Cit paper
        reports a similar number, 1266, for these three files).
        """
        self._prepare("cit", ["out"])
        expected = self._expected("cit")
        _run_rescoring("rescoring_cit.json")
        out_dir = DATA_DIR / "cit" / "out"
        fdr_dir = out_dir / "results" / "percolator"
        self._check_files(out_dir, expected["output_files"])
        self._check_columns(fdr_dir / "original.tab", expected["original_columns"])
        self._check_columns(fdr_dir / "rescore.tab", expected["rescore_columns"])
        loc_tab = fdr_dir / "localize_mod" / "rescore.tab"
        self._check_columns(loc_tab, expected["columns"]["localize_mod/rescore.tab"])
        # citrulline neutral loss (isocyanic acid) features
        self.assertIn("ANNOTATED_NL_COUNT", _tab_columns(loc_tab))

        psms = _accepted(fdr_dir / "localize_mod" / "rescore.percolator.psms.txt")
        self.assertTrue(psms["peptide"].str.contains(r"\[UNIMOD:7\]", regex=True).all())
        n_cit = psms["peptide"].str.contains(r"R\[UNIMOD:7\]", regex=True).sum()
        n_cit_reference = sum(expected["reference_citrullinated_psms_1pct_fdr"]["this_run_localize_mod"].values())
        self.assertGreater(n_cit, 0.8 * n_cit_reference, f"{n_cit} citrullinated PSMs, reference {n_cit_reference}")

    def test_rescoring_multifrag(self):
        """
        Rescoring of EID, UVPD and ETciD spectra with Prosit MultiFrag (PXD065289), one job per activation method.

        The activation of EID and UVPD spectra cannot be read from the mzML files (beam-type CID / "hcd" in the filter
        string), so every file is rescored with its own config that sets fragmentation_method and the ion types. ETciD
        is set as well, so the test does not depend on reading the supplemental activation from the mzML.
        """
        methods = {"eid": "aAbcxXyzZ", "uvpd": "aAbCyz", "etcid": "AbcCyzZ"}
        self._prepare("multifrag", [f"out_{method}" for method in methods])
        for method, ion_types in methods.items():
            with self.subTest(method=method):
                expected = self._expected("multifrag", f"expected_{method}.json")
                _run_rescoring(f"rescoring_multifrag_{method}.json")
                fdr_dir = self._check_standard_run(DATA_DIR / "multifrag" / f"out_{method}", expected)
                # MultiFrag features are calculated per configured ion type
                columns = _tab_columns(fdr_dir / "rescore.tab")
                for ion_type in ion_types:
                    self.assertIn(f"count_observed_and_predicted_{ion_type}", columns)

    def test_rescoring_phospho_localization(self):
        """
        Phospho site localization on a synthetic phosphopeptide library (PXD007058, pools 1-3, HCD OT).

        With ptm_localization, original.tab contains all localization candidates of every spectrum (including the
        decoy site candidates), so the original percolator results are no valid baseline. Instead, the localizations
        are checked against the known sites of the synthetic library.
        """
        self._prepare("phospho", ["out"])
        expected = self._expected("phospho")
        _run_rescoring("rescoring_phospho.json")
        out_dir = DATA_DIR / "phospho" / "out"
        fdr_dir = out_dir / "results" / "percolator"
        self._check_files(out_dir, expected["output_files"])
        self._check_columns(fdr_dir / "original.tab", expected["original_columns"])
        self._check_columns(fdr_dir / "rescore.tab", expected["rescore_columns"])
        loc_dir = fdr_dir / "localize_phospho"
        with_prob = loc_dir / "rescore.percolator_with_prob.psms.txt"
        self._check_columns(with_prob, expected["columns"]["localize_phospho/rescore.percolator_with_prob.psms.txt"])

        psms = _accepted(with_prob)
        self.assertGreater(len(psms), 0.5 * expected["reference_counts_1pct_fdr"][str(with_prob.relative_to(fdr_dir))])
        self.assertTrue(psms["peptide"].str.contains(r"[STY]\[UNIMOD:21\]", regex=True).all())
        self.assertTrue(psms["loc_prob"].between(0, 1).all())
        self.assertTrue(psms["localization_flr"].between(0, 1).all())

        # localizations at 1% estimated FLR agree with the known sites of the synthetic library
        truth = pd.read_csv(DATA_DIR / "phospho" / "expected_outputs" / "library_truth.tsv", sep="\t")
        true_sites = {
            (pool, seq): {int(p) for p in str(positions).split("&")}
            for pool, seq, positions in truth.itertuples(index=False)
        }
        confident = psms[psms["localization_flr"] <= FDR]
        correct = wrong = 0
        for filename, peptide in zip(confident["filename"], confident["peptide"], strict=True):
            sequence, sites = _phospho_sites(peptide)
            key = (re.search(r"(pool\d)", filename).group(1), sequence)
            if key in true_sites:
                correct += sites == true_sites[key]
                wrong += sites != true_sites[key]
        self.assertGreater(correct, 0)
        self.assertGreaterEqual(correct / (correct + wrong), 0.85, f"{wrong} of {correct + wrong} wrong localizations")

    def test_rescoring_hla_ptm_fdr(self):
        """
        Separate FDR control of modified PSMs on HLA class I peptides with the 22PTM models (PXD017448).

        With a PTM model and without localization, ptm_fdr_control is on by default: unmodified PSMs are rescored in
        results/percolator, modified PSMs (oxidation, phospho, deamidation, acetylation) separately in mod_fdr with
        0/1 indicator features per PTM.
        """
        self._prepare("hla", ["out"])
        expected = self._expected("hla")
        _run_rescoring("rescoring_hla.json")
        out_dir = DATA_DIR / "hla" / "out"
        fdr_dir = out_dir / "results" / "percolator"
        mod_dir = fdr_dir / "mod_fdr"
        self._check_files(out_dir, expected["output_files"])
        for tab in ("original.tab", "rescore.tab", "mod_fdr/original.tab", "mod_fdr/rescore.tab"):
            self._check_columns(fdr_dir / tab, expected["columns"][tab])
        self.assertTrue(any(column.startswith("ptm_UNIMOD_") for column in _tab_columns(mod_dir / "rescore.tab")))

        for level in ("psms", "peptides"):
            # unmodified and modified PSMs together
            self._check_gain(
                _n_accepted(fdr_dir / f"original.percolator.{level}.txt", mod_dir / f"original.percolator.{level}.txt"),
                _n_accepted(fdr_dir / f"rescore.percolator.{level}.txt", mod_dir / f"rescore.percolator.{level}.txt"),
                level,
            )
        self.assertFalse(_accepted(fdr_dir / "rescore.percolator.psms.txt")["peptide"].str.contains("UNIMOD").any())
        modified = _accepted(mod_dir / "rescore.percolator.psms.txt")["peptide"]
        self.assertGreater(len(modified), 0)
        self.assertTrue(modified.str.contains("UNIMOD").all())

    def test_rescoring_quantification(self):
        """
        MS1 quantification with FlashLFQ, picked group FDR and directLFQ on the ProteoBench Astral DDA ion module.

        Six Astral DDA runs (conditions A and B, three replicates each) of a human/yeast/E. coli mix, cut to a retention
        time window, searched with MSFragger DDA+ and rescored with chimeric peak subtraction. All quantification
        settings are measured from the data ("auto": MS1 tolerance, match-between-runs window, RT drift bin) and the
        intensities get the RT drift correction. Expected log2 A/B ratios: human 0, yeast 1, E. coli -2. FlashLFQ needs
        dotnet and its command line release (inputs/flashlfq/CMD.dll, see the integration test workflow); the test is
        skipped without them.
        """
        self._prepare("quant", ["out"])
        if not (DATA_DIR / "quant" / "inputs" / "flashlfq" / "CMD.dll").is_file() or shutil.which("dotnet") is None:
            self.skipTest("FlashLFQ (tests/data/quant/inputs/flashlfq/CMD.dll) or dotnet not installed.")
        expected = self._expected("quant")
        _run_rescoring("rescoring_quant.json")
        out_dir = DATA_DIR / "quant" / "out"
        quant_dir = out_dir / "quantification"
        fdr_dir = self._check_standard_run(out_dir, expected)
        self._check_files(quant_dir, expected["quantification_files"])
        ranges = expected["ranges"]

        # PSMs and peptides at 1% FDR before and after rescoring
        for level in ("psms", "peptides"):
            for scores, value_range in ranges[f"{level}_1pct_fdr"].items():
                n = _n_accepted(fdr_dir / f"{scores}.percolator.{level}.txt")
                self._check_range(n, value_range, f"{scores} {level} at 1% FDR")

        # every setting was measured from the data
        with open(quant_dir / "settings.json") as f:
            settings = json.load(f)
        for name, value_range in ranges["settings"].items():
            self.assertEqual(settings[name]["source"], "auto", name)
            self._check_range(settings[name]["value"], value_range, name)

        # FlashLFQ quantified peaks per run
        peaks = pd.read_csv(quant_dir / "flashlfq" / "QuantifiedPeaks.tsv", sep="\t", low_memory=False)
        quantified = peaks[peaks["Peak intensity"] > 0].groupby("File Name").size()
        for run, value_range in ranges["flashlfq_peaks"].items():
            self._check_range(int(quantified.get(run, 0)), value_range, f"FlashLFQ peaks in {run}")

        # picked group FDR protein groups at 1% FDR
        groups = pd.read_csv(quant_dir / "picked_group_fdr" / "proteinGroups.txt", sep="\t")
        n_groups = int(((groups["Q-value"] <= FDR) & groups["Reverse"].isna()).sum())
        self.assertGreater(n_groups, 0.8 * expected["reference"]["protein_groups_1pct_fdr"])

        # precursor (FlashLFQ) and protein (directLFQ) ratios
        ratios = {k: float(v) for k, v in expected["expected_log2_ratios"].items()}
        for table_name, reference in expected["reference"]["ratios"].items():
            table = pd.read_csv(quant_dir / table_name, sep="\t")
            result = _species_ratios(table, expected["condition_a"], expected["condition_b"], ratios)
            with self.subTest(table=table_name):
                self.assertGreater(result["entries"], 0.8 * reference["entries"], result)
                self.assertLess(result["median_abs_epsilon"], reference["median_abs_epsilon"] + 0.03, result)
                for name, ratio in ratios.items():
                    tolerance = expected["log2_ratio_tolerance"][name]
                    self.assertAlmostEqual(result[f"log2_{name}"], ratio, delta=tolerance, msg=f"{name}: {result}")
