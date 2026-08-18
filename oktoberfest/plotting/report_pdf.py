"""Render an investigation report (``investigate_*.html``) as a printable PDF.

Works on an ALREADY-BUILT report — no re-run and no percolator dir needed. The HTML is rewritten for
print and handed to headless Chrome:

  * the interactive spectra viewer is dropped (a Plotly widget with one visible trace and a dropdown;
    on paper it prints as a single arbitrary spectrum plus ~3 MB of dead JavaScript),
  * every collapsed ``<details>`` is expanded, so nothing is hidden in a static document,
  * print CSS is added: landscape pages, page breaks between sections, no sticky nav, figures kept
    off the page seams.

Everything is guarded: the PDF is a convenience on top of the HTML report, so a missing Chrome or a
failed render costs the PDF only. Run it from the command line for an existing report with
``python -m oktoberfest.plotting.report_pdf <report.html>``.
"""

import argparse
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

#: Chrome/Chromium executables tried, in order, for the headless render.
CHROME_CANDIDATES = ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome")

#: Wall-clock ceiling for the headless render. A bulk report is a large document, but a Chrome that
#: has not finished after this is hung, and a hung subprocess must not hold up a rescoring run.
RENDER_TIMEOUT_S = 900

PRINT_CSS = """
/* landscape: every non-native panel in this report is a wide multi-panel figure, and on portrait A4
   they shrink to about half the usable width with the rest of the page left blank. */
@page{size:A4 landscape;margin:10mm 10mm}
@media print{
  nav.toc{display:none!important}
  body{background:#fff}
  .wrap{max-width:none;padding:0}
  details.sec{break-inside:auto;margin:0 0 10px}
  details.sec>summary{break-after:avoid;break-inside:avoid}
  details.sec+details.sec>summary{break-before:page}
  .figcard,.svgcard,.kpi,.cfg,.masthead{break-inside:avoid}
  /* section header + figure + "what this shows" caption must fit one landscape page together,
     otherwise the caption orphans onto the next page away from the figure it explains */
  .figcard img{max-height:118mm;object-fit:contain}
  .svgcard img{max-height:105mm;object-fit:contain}
  .grid2{grid-template-columns:repeat(auto-fit,minmax(115mm,1fr))}
  .figcap{font-size:10.5px;line-height:1.45}
  details.raw{break-inside:avoid}
  footer{break-before:avoid}
}
details.sec>summary::after{content:''!important}
"""


def find_chrome() -> Optional[str]:
    """Path to a usable headless Chrome/Chromium, or None if none is installed."""
    for candidate in CHROME_CANDIDATES:
        found = shutil.which(candidate)
        if found:
            return found
    return None


def drop_section(html: str, sid: str):
    """Remove one ``<details class='sec' id='sid'>…</details>`` block by tag counting.

    The report nests ``<details>`` (the per-raw gallery holds one per raw file), so a regex for the
    closing tag would stop at the first inner ``</details>``; the nesting depth is counted instead.

    :param html: the report document
    :param sid: id of the section to remove
    :return: (document, whether the section was found and removed)
    """
    m = re.search(rf"<details class='sec'[^>]*id='{sid}'[^>]*>", html)
    if not m:
        return html, False
    i, depth = m.end(), 1
    tag = re.compile(r"<details\b|</details>")
    while depth:
        t = tag.search(html, i)
        if not t:
            return html, False  # unbalanced: leave the document alone rather than truncate it
        depth += 1 if t.group() == "<details" else -1
        i = t.end()
    return html[: m.start()] + html[i:], True


def prepare(html: str, keep_gallery: bool = True):
    """Rewrite a report for print.

    :param html: the report document
    :param keep_gallery: keep the per-raw-file appendix (3 SVGs per raw file, i.e. most of the pages)
    :return: (print-ready document, list of notes on what was changed)
    """
    notes = []
    html, ok = drop_section(html, "spectra")
    notes.append("spectra viewer dropped" if ok else "no spectra section found")
    if ok:  # the footer explains a viewer this document no longer contains
        html = html.replace(
            "Spectra are observed mzML vs a clean Koina re-query (not the rescoring's stored predictions). ", ""
        )
    if not keep_gallery:
        html, ok = drop_section(html, "gallery")
        notes.append("per-raw gallery dropped" if ok else "no gallery section found")
    # the viewer is gone, so strip the Plotly bundle it left behind (~3 MB of inline JS)
    html, n = re.subn(r"<script[^>]*>\s*/\*\*\s*plotly\.js.*?</script>", "", html, flags=re.S | re.I)
    if n:
        notes.append(f"stripped {n} plotly bundle(s)")
    # expand everything: a printed document has no affordance for clicking a disclosure triangle
    html = re.sub(r"<details(?![^>]*\bopen\b)", "<details open", html)
    # nav links are dead on paper, and the hash-opener script has nothing left to open
    html = re.sub(r"<nav class='toc'>.*?</nav>", "", html, flags=re.S)
    html = re.sub(r"<script>function _openHash.*?</script>", "", html, flags=re.S)
    html = html.replace("</style>", PRINT_CSS + "</style>", 1)
    return html, notes


def html_to_pdf(
    report: Union[str, Path], out_pdf: Optional[Union[str, Path]] = None, keep_gallery: bool = True, log=None
) -> Optional[Path]:
    """Render a built HTML report to PDF with headless Chrome.

    :param report: path of the HTML report
    :param out_pdf: where to write the PDF; defaults to the report path with a .pdf suffix
    :param keep_gallery: keep the per-raw-file appendix
    :param log: where progress goes; defaults to this module's logger
    :return: the path of the written PDF, or None if it could not be rendered
    """
    log = log or logger.info
    report = Path(report).resolve()
    out = Path(out_pdf).resolve() if out_pdf else report.with_suffix(".pdf")
    tmpdir = None
    try:
        chrome = find_chrome()
        if chrome is None:
            log(f"[pdf] no headless Chrome found (tried: {', '.join(CHROME_CANDIDATES)}), skipping the PDF")
            return None
        html, notes = prepare(report.read_text(), keep_gallery=keep_gallery)
        tmpdir = Path(tempfile.mkdtemp(prefix="okt_pdf_"))
        tmp = tmpdir / "print.html"
        tmp.write_text(html)
        log(f"[pdf] {report.name}: {', '.join(notes)} "
            f"({report.stat().st_size / 1e6:.0f} MB -> {len(html) / 1e6:.0f} MB)")
        cmd = [chrome, "--headless", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
               "--virtual-time-budget=120000", f"--print-to-pdf={out}", tmp.as_uri()]
        # noqa justified: every element of cmd is either a constant of this module or a path we wrote
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=RENDER_TIMEOUT_S)  # noqa: S603
        if not out.exists():
            log(f"[pdf] chrome failed to write the PDF, skipping it:\n{r.stderr[-2000:]}")
            return None
        log(f"[pdf] wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
        return out
    except Exception as e:  # noqa: BLE001 - the PDF is a convenience; it may never break a run
        log(f"[pdf] PDF rendering failed, skipped: {e}")
        return None
    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    """Render an existing report to PDF from the command line."""
    ap = argparse.ArgumentParser(description="Render an Oktoberfest investigation report as a printable PDF.")
    ap.add_argument("report", help="investigate_*.html built by oktoberfest.plotting.investigate")
    ap.add_argument("out_pdf", nargs="?", default=None)
    ap.add_argument("--no-gallery", action="store_true", help="also drop the per-raw-file appendix")
    a = ap.parse_args()
    if html_to_pdf(a.report, a.out_pdf, keep_gallery=not a.no_gallery, log=print) is None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
