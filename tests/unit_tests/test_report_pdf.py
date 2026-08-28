import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from oktoberfest.plotting import report_pdf as rp

SECTION = "<details class='sec' id='{sid}'><summary>{sid}</summary>{body}</details>"


def report(body_sections: str, style: str = "body{color:#000}") -> str:
    """Wrap sections in the minimal report chrome the print rewriter expects."""
    return (
        "<!doctype html><html><head><style>"
        + style
        + "</style></head><body><nav class='toc'><a href='#sa'>SA</a></nav>"
        + body_sections
        + "<script>function _openHash(){/* jump to a section */}</script></body></html>"
    )


class TestFindChrome(unittest.TestCase):
    """Locating a headless Chrome, and coping with its absence."""

    def test_returns_none_when_nothing_is_installed(self):
        """No Chrome on the PATH is a normal outcome, not an error."""
        with patch.object(rp.shutil, "which", return_value=None):
            self.assertIsNone(rp.find_chrome())

    def test_returns_the_first_candidate_found(self):
        """Candidates are tried in order, so the preferred binary wins."""
        with patch.object(rp.shutil, "which", side_effect=lambda c: "/usr/bin/chromium" if c == "chromium" else None):
            self.assertEqual(rp.find_chrome(), "/usr/bin/chromium")


class TestDropSection(unittest.TestCase):
    """Removing one section by counting <details> depth."""

    def test_removes_the_named_section_only(self):
        """The requested section goes; its neighbours stay untouched."""
        html = SECTION.format(sid="sa", body="keep") + SECTION.format(sid="spectra", body="drop")
        out, ok = rp.drop_section(html, "spectra")
        self.assertTrue(ok)
        self.assertNotIn("drop", out)
        self.assertIn("keep", out)

    def test_absent_section_is_a_no_op(self):
        """A report built without the spectra viewer must survive unchanged."""
        html = SECTION.format(sid="sa", body="keep")
        out, ok = rp.drop_section(html, "spectra")
        self.assertFalse(ok)
        self.assertEqual(out, html)

    def test_nested_details_do_not_truncate_the_document(self):
        """The gallery nests one <details> per raw file; a naive regex would stop at the first one."""
        nested = "<details><summary>raw_00</summary>inner</details>"
        html = (
            SECTION.format(sid="gallery", body=nested * 3)
            + "<p id='after'>after</p>"
            + SECTION.format(sid="sa", body="keep")
        )
        out, ok = rp.drop_section(html, "gallery")
        self.assertTrue(ok)
        self.assertNotIn("inner", out)
        self.assertIn("after", out)
        self.assertIn("keep", out)

    def test_unbalanced_markup_leaves_the_document_alone(self):
        """Truncating a report is worse than not shrinking it."""
        html = "<details class='sec' id='spectra'>never closed<p>tail</p>"
        out, ok = rp.drop_section(html, "spectra")
        self.assertFalse(ok)
        self.assertEqual(out, html)


class TestPrepare(unittest.TestCase):
    """The HTML-to-print rewrite, which runs with or without Chrome present."""

    def _prepared(self, **kwargs):
        html = report(
            SECTION.format(sid="spectra", body="<script>/** plotly.js v2 */ var x=1;</script>viewer")
            + SECTION.format(sid="gallery", body="<details><summary>raw</summary>PER_RAW_PANEL</details>")
            + SECTION.format(sid="sa", body="panel")
        )
        return rp.prepare(html, **kwargs)

    def test_spectra_viewer_and_its_plotly_bundle_are_dropped(self):
        """The viewer prints as one arbitrary spectrum plus megabytes of dead JS."""
        out, notes = self._prepared()
        self.assertNotIn("viewer", out)
        self.assertNotIn("plotly.js", out)
        self.assertIn("spectra viewer dropped", notes)

    def test_gallery_is_kept_by_default_and_droppable(self):
        """The per-raw appendix is most of the pages, so dropping it is opt-in."""
        self.assertIn("PER_RAW_PANEL", self._prepared()[0])
        self.assertNotIn("PER_RAW_PANEL", self._prepared(keep_gallery=False)[0])

    def test_every_details_is_expanded(self):
        """Nothing may stay hidden behind a disclosure triangle in a static document."""
        out, _ = self._prepared()
        self.assertNotIn("<details class=", out.replace("<details open class=", ""))
        self.assertEqual(out.count("<details"), out.count("<details open"))

    def test_dead_navigation_is_removed_and_print_css_added(self):
        """Nav links and the hash-opener do nothing on paper; the print CSS must land inside <style>."""
        out, _ = self._prepared()
        self.assertNotIn("nav class='toc'", out)
        self.assertNotIn("_openHash", out)
        self.assertIn("@page", out)
        self.assertLess(out.index("@page"), out.index("</style>"))

    def test_a_report_without_a_viewer_is_still_print_ready(self):
        """Reports built with "spectra": false must not be treated as malformed."""
        out, notes = rp.prepare(report(SECTION.format(sid="sa", body="panel")))
        self.assertIn("no spectra section found", notes)
        self.assertIn("@page", out)
        self.assertIn("panel", out)


class TestHtmlToPdf(unittest.TestCase):
    """The PDF is a convenience on top of the HTML; nothing about it may break a run."""

    @staticmethod
    def _report(tmp: str) -> Path:
        path = Path(tmp) / "investigate_report.html"
        path.write_text(report(SECTION.format(sid="sa", body="panel")))
        return path

    def test_missing_chrome_returns_none_and_says_so(self):
        """No Chrome is a skip with an explanation, not a crash."""
        messages = []
        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value=None):
            self.assertIsNone(rp.html_to_pdf(self._report(tmp), log=messages.append))
        self.assertTrue(any("no headless Chrome" in m for m in messages))

    def test_chrome_that_writes_nothing_returns_none(self):
        """A silent Chrome failure must not leave the caller believing a PDF exists."""
        messages = []
        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value="/usr/bin/chrome"):
            with patch.object(rp.subprocess, "run") as run:
                run.return_value = type("R", (), {"stderr": "boom", "stdout": "", "returncode": 1})()
                self.assertIsNone(rp.html_to_pdf(self._report(tmp), log=messages.append))
        self.assertTrue(any("chrome failed" in m for m in messages))

    def test_a_raising_chrome_is_swallowed(self):
        """Any exception below this call costs the PDF only."""
        messages = []
        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value="/usr/bin/chrome"):
            with patch.object(rp.subprocess, "run", side_effect=OSError("no exec")):
                self.assertIsNone(rp.html_to_pdf(self._report(tmp), log=messages.append))
        self.assertTrue(any("PDF rendering failed" in m for m in messages))

    def test_written_pdf_is_returned_and_the_tempdir_cleaned_up(self):
        """On success the caller gets the path, and no scratch directory is left behind."""
        seen = {}

        def fake_run(cmd, **kwargs):
            seen["cmd"] = cmd
            out = next(a.split("=", 1)[1] for a in cmd if a.startswith("--print-to-pdf="))
            Path(out).write_bytes(b"%PDF-1.4\n")
            seen["tmp_html"] = Path(cmd[-1].replace("file://", ""))
            return type("R", (), {"stderr": "", "stdout": "", "returncode": 0})()

        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value="/usr/bin/chrome"):
            with patch.object(rp.subprocess, "run", side_effect=fake_run):
                result = rp.html_to_pdf(self._report(tmp), log=lambda _: None)
            self.assertEqual(result, Path(tmp) / "investigate_report.pdf")
            self.assertTrue(result.exists())
        self.assertIn("--headless", seen["cmd"])
        self.assertFalse(seen["tmp_html"].parent.exists(), "scratch directory was not cleaned up")

    def test_a_stale_pdf_is_not_reported_as_freshly_written(self):
        """Chrome failing on a re-render must not hand back the previous run's PDF."""
        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value="/usr/bin/chrome"):
            report = self._report(tmp)
            stale = report.with_suffix(".pdf")
            stale.write_bytes(b"%PDF-stale")

            def failing_run(cmd, **kwargs):
                return type("R", (), {"stderr": "boom", "stdout": "", "returncode": 1})()

            with patch.object(rp.subprocess, "run", side_effect=failing_run):
                self.assertIsNone(rp.html_to_pdf(report, log=lambda _: None))
            self.assertFalse(stale.exists(), "the stale PDF should have been cleared, not returned")

    def test_out_path_override_is_honoured(self):
        """A caller-supplied destination must win over the default .pdf sibling."""
        with TemporaryDirectory() as tmp, patch.object(rp, "find_chrome", return_value="/usr/bin/chrome"):
            target = Path(tmp) / "elsewhere.pdf"

            def fake_run(cmd, **kwargs):
                Path(next(a.split("=", 1)[1] for a in cmd if a.startswith("--print-to-pdf="))).write_bytes(b"%PDF")
                return type("R", (), {"stderr": "", "stdout": "", "returncode": 0})()

            with patch.object(rp.subprocess, "run", side_effect=fake_run):
                self.assertEqual(rp.html_to_pdf(self._report(tmp), target, log=lambda _: None), target)
