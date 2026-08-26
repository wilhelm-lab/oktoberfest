"""Tests for the HTML/PDF investigation report and the config option that switches it on."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import matplotlib
import numpy as np

matplotlib.use("Agg")

import oktoberfest.plotting as pl  # noqa: E402
from oktoberfest.plotting import investigate as inv  # noqa: E402
from oktoberfest.utils import Config  # noqa: E402
from oktoberfest.utils.config import REPORT_DEFAULTS  # noqa: E402

TAB_COLUMNS = [
    "SpecId", "Label", "ScanNr", "filename", "ExpMass", "Mass", "RT", "abs_rt_diff",
    "collision_energy_aligned", "missedCleavages", "sequence_length", "KR", "delta_mass_ppm",
    "pearson_corr", "spectral_angle", "fraction_observed_and_predicted", "Peptide",
    "Charge1", "Charge2", "Charge3",
]


def write_run(root: Path, n: int = 400, n_files: int = 6, method: str = "percolator", ragged: bool = True) -> Path:
    """Write a minimal rescoring output tree the report can be built from.

    :param root: directory to write the run into
    :param n: number of PSMs
    :param n_files: number of raw files the PSMs are spread over
    :param method: "percolator" or "mokapot" — file names and column names follow it
    :param ragged: give the output files a trailing protein column with extra tabs, as percolator does
    :return: the run's out/ directory
    """
    rng = np.random.default_rng(0)
    fdr_dir = root / "out" / "results" / method
    fdr_dir.mkdir(parents=True)
    label = np.where(rng.random(n) < 0.4, -1, 1)
    charge = rng.integers(2, 4, n)
    spectral_angle = np.clip(rng.beta(np.where(label == 1, 5, 2), 3, n), 0, 1)
    score = rng.normal(np.where(label == 1, 1.5, -1.0), 0.8, n)
    peptides = np.array([f"PEPTIDEK{i % 40}" for i in range(n)])
    files = [f"raw_{i % n_files:02d}" for i in range(n)]
    spec_ids = [f"{files[i]}-{i}-{peptides[i]}-{charge[i]}" for i in range(n)]

    with open(fdr_dir / "rescore.tab", "w") as f:
        f.write("\t".join([*TAB_COLUMNS, "Proteins"]) + "\n")
        for i in range(n):
            row = [spec_ids[i], label[i], (i + 1) * 2, files[i], i, 800 + i, i % 60, abs(rng.normal(2, 1)),
                   0.3, i % 3, len(peptides[i]), 1, rng.normal(0, 3), rng.random(), spectral_angle[i],
                   rng.random(), f"_.{peptides[i]}._", *(int(charge[i] == c) for c in (1, 2, 3))]
            proteins = "\t".join(f"sp|P{i:04d}|PROT{j}" for j in range(1 + (i % 3 if ragged else 0)))
            f.write("\t".join(str(v) for v in row) + "\t" + proteins + "\n")

    id_col, score_col, q_col, pep_col = (
        ("PSMId", "score", "q-value", "peptide") if method == "percolator"
        else ("SpecId", "mokapot score", "mokapot q-value", "Peptide")
    )

    def qvalues(scores):
        """Textbook target-decoy q-values, so the report's counts are checkable by hand."""
        order = np.argsort(-scores)
        fdr = np.cumsum(label[order] == -1) / np.maximum(np.cumsum(label[order] == 1), 1)
        q = np.empty(n)
        q[order] = np.minimum.accumulate(np.clip(fdr, 0, 1)[::-1])[::-1]
        return q

    for prefix, scores in (("rescore", score), ("original", score - 0.3)):
        q = qvalues(scores)
        for lab, name in ((1, "psms"), (-1, "decoy.psms"), (1, "peptides"), (-1, "decoy.peptides")):
            rows = np.where(label == lab)[0]
            rows = rows[np.argsort(-scores[rows])]
            if "peptides" in name:  # one row per distinct peptide, best first
                _, first = np.unique(peptides[rows], return_index=True)
                rows = rows[np.sort(first)]
            with open(fdr_dir / f"{prefix}.{method}.{name}.txt", "w") as f:
                f.write(f"{id_col}\t{score_col}\t{q_col}\tposterior_error_prob\t{pep_col}\tproteinIds\n")
                for i in rows:
                    trailing = "\tsp|Q{i:04d}|EXTRA" if ragged else ""
                    f.write(f"{spec_ids[i]}\t{scores[i]:.6f}\t{q[i]:.6g}\t0.1\t_.{peptides[i]}._\t"
                            f"sp|P{i:04d}|PROT{trailing}\n")

    svg = "<svg xmlns='http://www.w3.org/2000/svg' width='80' height='40'><text x='4' y='20'>{}</text></svg>"
    for stem in ("psm_1%_FDR", "peptide_1%_FDR", "target_vs_decoys_sa_distribution"):
        (fdr_dir / f"{stem}.svg").write_text(svg.format(stem))
    return root / "out"


class TestReportConfig(unittest.TestCase):
    """The config option that switches the report on."""

    @staticmethod
    def _config(data):
        config = Config()
        config.data = data
        return config

    def test_disabled_by_default(self):
        """A config that does not mention the report does not get one."""
        self.assertEqual(self._config({}).report, REPORT_DEFAULTS)
        self.assertFalse(self._config({}).report["enabled"])

    def test_boolean_shorthand(self):
        """"report": true enables the report and leaves every other setting at its default."""
        settings = self._config({"report": True}).report
        self.assertTrue(settings["enabled"])
        self.assertEqual(settings["n_per_group"], REPORT_DEFAULTS["n_per_group"])

    def test_dictionary_overrides(self):
        """Individual settings can be overridden, and are coerced to the type of their default."""
        settings = self._config({"report": {"enabled": 1, "spectra": False, "n_per_group": "5"}}).report
        self.assertIs(settings["enabled"], True)
        self.assertIs(settings["spectra"], False)
        self.assertEqual(settings["n_per_group"], 5)

    def test_rejects_malformed_options(self):
        """A malformed option is refused loudly: silently not writing a report is the worse failure."""
        for data in ({"report": "yes"}, {"report": {"enable": True}}, {"report": {"n_per_group": "many"}}):
            with self.assertRaises(ValueError):
                _ = self._config(data).report

    def test_checked_before_the_run(self):
        """Config.check() validates the option, so a typo fails before the compute, not after it."""
        config = self._config({"report": {"enable": True}, "models": {"intensity": "x", "irt": "y"}, "tag": ""})
        with patch.object(Config, "_check_tmt"), patch.object(Config, "_check_koina_model_availability"), \
                patch.object(Config, "job_type", "Rescoring"), patch.object(Config, "quantification", False):
            with self.assertRaises(ValueError):
                config.check()


class TestReportGate(unittest.TestCase):
    """The wiring between the config option and report generation."""

    def test_not_written_unless_enabled(self):
        """No report is built when the option is absent or off; it is built when it is on."""
        with TemporaryDirectory() as tmp:
            data_dir = write_run(Path(tmp)) / "results" / "percolator"
            config = Config()
            for data, expected in (({}, False), ({"report": False}, False), ({"report": {"enabled": True}}, True)):
                config.data = data
                with patch("oktoberfest.plotting.investigate.build_report_safe") as build:
                    pl.plot_investigation_report(data_dir, config)
                self.assertEqual(build.called, expected, f"config {data}")

    def test_settings_are_forwarded(self):
        """Every report setting reaches the builder — an option that is not wired through is a no-op."""
        with TemporaryDirectory() as tmp:
            data_dir = write_run(Path(tmp)) / "results" / "percolator"
            config = Config()
            config.data = {"report": {"enabled": True, "spectra": False, "n_per_group": 3, "pdf": True,
                                      "max_psms": 10, "max_gallery_files": 2, "max_embedded_mb": 5}}
            with patch("oktoberfest.plotting.investigate.build_report_safe") as build:
                pl.plot_investigation_report(data_dir, config)
            kwargs = build.call_args.kwargs
            self.assertEqual(kwargs["want_spectra"], False)
            self.assertEqual(kwargs["n_per_group"], 3)
            self.assertEqual(kwargs["want_pdf"], True)
            self.assertEqual(kwargs["max_psms"], 10)
            self.assertEqual(kwargs["max_gallery_files"], 2)
            self.assertEqual(kwargs["max_embedded_mb"], 5)

    def test_a_broken_report_never_breaks_the_run(self):
        """The report is a diagnostic aid: whatever it does, the rescoring run must survive it."""
        config = Config()
        config.data = {"report": True}
        with patch("oktoberfest.plotting.investigate.build_report_safe", side_effect=RuntimeError("boom")):
            self.assertIsNone(pl.plot_investigation_report(Path("/nonexistent"), config))
        config.data = {"report": "malformed"}
        self.assertIsNone(pl.plot_investigation_report(Path("/nonexistent"), config))


class TestReportBuild(unittest.TestCase):
    """Building the report itself."""

    def test_percolator_run(self):
        """A percolator run yields a self-contained HTML with the sections its data supports."""
        with TemporaryDirectory() as tmp:
            out_dir = write_run(Path(tmp))
            report = inv.build_report(out_dir, want_spectra=False, log=lambda m: None)
            self.assertTrue(report.exists())
            html = report.read_text()
            for section in ("summary", "yield", "movement", "sa", "headroom", "calibration", "native"):
                self.assertIn(f"id='{section}'", html, f"missing section {section}")
            self.assertNotIn("http://", html.replace("http://www.w3.org", ""))  # self-contained

    def test_mokapot_run(self):
        """A mokapot run is read too: its files and column names differ from percolator's."""
        with TemporaryDirectory() as tmp:
            out_dir = write_run(Path(tmp), method="mokapot")
            paths = inv.find_paths(out_dir)
            self.assertEqual(paths["method"], "mokapot")
            report = inv.build_report(out_dir, want_spectra=False, log=lambda m: None)
            self.assertIn("id='summary'", report.read_text())

    def test_counts_match_the_fdr_output(self):
        """The headline numbers are the accepted rows of the FDR output, counted independently here."""
        with TemporaryDirectory() as tmp:
            out_dir = write_run(Path(tmp))
            paths = inv.find_paths(out_dir)
            stats = inv.overview_stats(paths, inv.load_data(paths))

            def accepted(name, key):
                with open(paths["fdr_dir"] / name) as f:
                    header = f.readline().rstrip("\n").split("\t")
                    ki, qi = header.index(key), header.index("q-value")
                    return {p[ki] for p in (line.rstrip("\n").split("\t") for line in f) if float(p[qi]) <= 0.01}

            rescored = accepted("rescore.percolator.psms.txt", "PSMId")
            original = accepted("original.percolator.psms.txt", "PSMId")
            self.assertEqual(stats["acc_psm"], len(rescored))
            self.assertEqual(stats["acc_pep"], len(accepted("rescore.percolator.peptides.txt", "peptide")))
            self.assertEqual(stats["gained_psm"], len(rescored - original))
            self.assertEqual(stats["lost_psm"], len(original - rescored))

    def test_missing_optional_columns(self):
        """Feature columns the run did not produce cost their own panels, not the report."""
        with TemporaryDirectory() as tmp:
            out_dir = write_run(Path(tmp))
            tab = out_dir / "results" / "percolator" / "rescore.tab"
            lines = tab.read_text().splitlines()
            drop = lines[0].split("\t").index("spectral_angle")
            tab.write_text("\n".join("\t".join(c for i, c in enumerate(ln.split("\t")) if i != drop)
                                     for ln in lines) + "\n")
            html = inv.build_report(out_dir, want_spectra=False, log=lambda m: None).read_text()
            self.assertIn("id='summary'", html)
            self.assertNotIn("id='sa'", html)  # the SA panels have nothing to show

    def test_guard_rails(self):
        """The size guard rails cap what is parsed and what is embedded."""
        with TemporaryDirectory() as tmp:
            out_dir = write_run(Path(tmp), n=400, n_files=6)
            with self.assertRaises(inv.ReportTooLargeError):
                inv.build_report(out_dir, want_spectra=False, max_psms=10, log=lambda m: None)
            self.assertIsNone(inv.build_report_safe(out_dir, want_spectra=False, max_psms=10, log=lambda m: None))
            results = out_dir / "results"
            for raw in {f"raw_{i:02d}" for i in range(6)}:
                (results / f"{raw}_irt_vs_pred_rt.svg").write_text("<svg xmlns='http://www.w3.org/2000/svg'/>")
            _, gallery, note = inv.collect_svgs(inv.find_paths(out_dir), max_gallery_files=2, max_embedded_mb=60)
            self.assertEqual(len(gallery), 2)
            self.assertIn("2 of 6 raw files", note)

    def test_both_parsers_agree(self):
        """The fast parser and the fallback must produce the same arrays, or results depend on pandas."""
        with TemporaryDirectory() as tmp:
            tab = write_run(Path(tmp)) / "results" / "percolator" / "rescore.tab"
            with open(tab) as f:
                header = f.readline().rstrip("\n").split("\t")
            cols = ["SpecId", "Label", "Mass", "spectral_angle", "filename", "Peptide"]
            idx = [header.index(c) for c in cols]
            parsed = []
            for chunks in (inv._chunks_pandas(tab, header, cols, idx), inv._chunks_python(tab, cols, idx)):
                maps = {c: {} for c in cols if c in inv.CAT_COLS}
                parts = {c: [] for c in cols}
                for raw in chunks:
                    for c in cols:
                        parts[c].append(inv._column_chunk(c, raw[c], maps.get(c)))
                parsed.append({c: np.concatenate(parts[c]) for c in cols})
            for c in cols:
                self.assertTrue(np.array_equal(parsed[0][c], parsed[1][c]), f"parsers disagree on {c}")


class TestRunParameters(unittest.TestCase):
    """Reading the experiment's settings off the run's config instead of assuming them.

    The report has to be as correct for a label-free bulk run as for a labelled single-cell one, so
    nothing that is a property of the SAMPLE may be a constant in the code.
    """

    #: Monoisotopic reporter-ion m/z that each label's masking window has to cover.
    REPORTERS = {"tmt": (126.1277, 131.1445), "tmtpro": (126.1277, 134.1548),
                 "itraq4": (114.1112, 117.1150), "itraq8": (113.1078, 121.1220)}

    def test_a_label_free_run_masks_nothing(self):
        """Without a tag there is no reporter region: that m/z range holds real fragment ions."""
        params = inv.run_params({"type": "Rescoring", "models": {"intensity": "Prosit_2020_intensity_HCD"}})
        self.assertIsNone(params["reporter"])
        self.assertEqual(params["tag"], "")
        self.assertIn("label-free", params["label"])

    def test_each_label_masks_its_own_reporter_region(self):
        """Every supported tag gets a window that brackets that plex's reporter ions."""
        for tag, (lowest, highest) in self.REPORTERS.items():
            for spelling in (tag, tag.upper(), f"{tag}_msa"):  # the tags Oktoberfest accepts
                lo, hi = inv.run_params({"tag": spelling})["reporter"]
                self.assertLess(lo, lowest, f"{spelling} window starts above its first reporter")
                self.assertGreater(hi, highest, f"{spelling} window ends below its last reporter")

    def test_a_label_masks_only_its_own_region(self):
        """The window is not a blanket low-m/z cut: b2/y1-sized fragments below it survive."""
        lo, _ = inv.run_params({"tag": "tmt"})["reporter"]
        self.assertGreater(lo, 121.2, "the TMT window reaches into the iTRAQ/immonium region")

    def test_tolerance_is_the_one_the_run_annotated_with(self):
        """An explicit massTolerance wins; otherwise the mass-analyzer fallback the run itself used."""
        explicit = inv.run_params({"massTolerance": 10, "unitMassTolerance": "ppm"})
        self.assertEqual((explicit["tol"], explicit["tol_unit"]), (10.0, "ppm"))
        self.assertIn("massTolerance", explicit["tol_source"])
        for analyzer, expected in (("FTMS", (20.0, "ppm")), ("TOF", (40.0, "ppm")), ("ITMS", (0.35, "da"))):
            params = inv.run_params({"inputs": {"instrument_type": analyzer}})
            self.assertEqual((params["tol"], params["tol_unit"]), expected, analyzer)
            self.assertIn(analyzer, params["tol_source"])
        unknown = inv.run_params({})  # no analyzer stated anywhere: FTMS, and said to be an assumption
        self.assertEqual((unknown["tol"], unknown["tol_unit"]), (20.0, "ppm"))
        self.assertIn("assumed", unknown["tol_source"])

    def test_a_half_specified_tolerance_falls_back(self):
        """A massTolerance without its unit is not usable — the same rule the annotation follows."""
        params = inv.run_params({"massTolerance": 10, "inputs": {"instrument_type": "TOF"}})
        self.assertEqual((params["tol"], params["tol_unit"]), (40.0, "ppm"))

    def test_match_window_follows_the_unit(self):
        """A ppm window scales with m/z, a dalton window does not."""
        predicted = np.array([200.0, 200.004, 200.02])
        ppm = inv.run_params({"massTolerance": 20, "unitMassTolerance": "ppm"})
        dalton = inv.run_params({"massTolerance": 0.01, "unitMassTolerance": "da"})
        self.assertEqual(inv.match_indices(200.0, predicted, ppm).tolist(), [0, 1])  # 20 ppm = 0.004 Da
        self.assertEqual(inv.match_indices(200.0, predicted, dalton).tolist(), [0, 1])
        self.assertEqual(inv.match_indices(2000.0, np.array([2000.0, 2000.02]), ppm).tolist(), [0, 1])
        self.assertEqual(inv.match_indices(2000.0, np.array([2000.0, 2000.02]), dalton).tolist(), [0])

    def test_koina_requery_uses_the_run_settings(self):
        """The spectra viewer re-queries with the run's own fragmentation method and server."""
        params = inv.run_params({"fragmentation_method": "hcd", "prediction_server": "example.org:443",
                                 "ssl": False})
        self.assertEqual(params["fragmentation"], "HCD")
        self.assertEqual(params["server"], "example.org:443")
        self.assertIs(params["ssl"], False)

    def test_parameter_table_mirrors_the_config(self):
        """Rows are keyed by their config path, and a value the config did not set is marked."""
        cfg = {"type": "Rescoring", "tag": "tmt", "massTolerance": 20, "unitMassTolerance": "ppm",
               "inputs": {"search_results_type": "msfragger"}, "matching_method": "global_ransac"}
        rows = {key: (value, note) for _, _, key, value, note in inv.config_rows(cfg, inv.run_params(cfg))}
        self.assertEqual(rows["inputs.search_results_type"], ("msfragger", ""))
        self.assertEqual(rows["matching_method"], ("global_ransac", ""))
        self.assertEqual(rows["massTolerance / unitMassTolerance"], ("20 ppm", ""))
        self.assertEqual(rows["fdr_estimation_method"][1], "default")  # not set -> marked, not invented
        self.assertEqual(rows["fdr_estimation_method"][0], "percolator")

    def test_parameter_table_drops_what_the_run_did_not_do(self):
        """Fasta digestion and custom modification rows appear only for a run that has them."""
        plain = {"type": "Rescoring"}
        keys = {key for _, _, key, _, _ in inv.config_rows(plain, inv.run_params(plain))}
        self.assertFalse(any(k.startswith("fastaDigestOptions") for k in keys))
        self.assertFalse(any("custom_modifications" in k for k in keys))
        digested = {"type": "Rescoring", "fastaDigestOptions": {"enzyme": "lysc"}}
        rows = {key: (value, note) for _, _, key, value, note in inv.config_rows(digested, inv.run_params(digested))}
        self.assertEqual(rows["fastaDigestOptions.enzyme"], ("lysc", ""))
        self.assertEqual(rows["fastaDigestOptions.missedCleavages"], ("2", "default"))

    def test_label_free_report_says_label_free(self):
        """End to end: a bulk run's report states the absence of a label rather than implying one."""
        cfg = {"type": "Rescoring", "inputs": {"search_results_type": "msfragger", "instrument_type": "TOF"},
               "models": {"intensity": "Prosit_2020_intensity_HCD"}}
        table = inv.config_table(cfg, inv.run_params(cfg))
        self.assertIn("label-free", table)
        self.assertNotIn("TMT", table)
        self.assertIn("40 ppm", table)  # the TOF tolerance the run annotated with
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = write_run(root)
            (root / "config.json").write_text(json.dumps(cfg))
            html = inv.build_report(out_dir, want_spectra=False, log=lambda m: None).read_text()
            self.assertIn("label-free", html)

    def test_report_without_a_config_admits_it(self):
        """No config.json means the parameters are unknown — not that they are the defaults."""
        table = inv.config_table({}, inv.run_params({}))
        self.assertIn("No config.json", table)
        self.assertNotIn("<table>", table)


if __name__ == "__main__":
    unittest.main()
