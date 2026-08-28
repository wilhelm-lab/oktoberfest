#!/usr/bin/env python3
"""Build the self-contained HTML investigation report for a finished rescoring run."""

# Point this at a Rescoring output folder (<run>/out, or the percolator/mokapot dir directly) and it
# writes one HTML file: identification yield vs FDR with and without Prosit, per-PSM score movement,
# identifications per raw file, percolator feature weights, spectral-angle diagnostics, rescoring
# headroom, RT and mass-error calibration, the native Oktoberfest plots plus a per-raw SVG gallery,
# and an interactive observed-vs-predicted spectra viewer.
#
# Generic by construction: percolator and mokapot output are both read, every optional feature column
# and SVG is used only if present, and the isobaric label, mass tolerance, fragmentation method and
# prediction server all come from the run's config.json rather than being assumed.
#
# Written during a run when the "report" config option is on (oktoberfest.utils.config.REPORT_DEFAULTS
# holds the guard rails). The caps only thin what is DRAWN or EMBEDDED -- every count, median and
# yield number is computed on all of the data.
#
# Two rules the panels are held to: a legend or annotation defines every colour, mark, unit and n; and
# a panel that restates a native plot, a KPI card or its neighbour does not go in (deliberate
# omissions are recorded where they would otherwise be re-added, see fig_sa_features/fig_calibration).
#
# Also runs standalone, auto-skipping the spectra viewer when mzML or Koina is unavailable:
#     python -m oktoberfest.plotting.investigate <run>/out [out.html] [--no-spectra] [--pdf]

import argparse
import base64
import contextlib
import io
import json
import re
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=False)  # force=False: never steal a backend the caller already set up
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

# ============================== house style ==============================
# The report chrome (masthead etc.) carries the Oktoberfest look; PLOT SERIES colours prioritise being
# easy to tell apart (a distinguishable qualitative palette), and the gains plot matches the NATIVE
# Oktoberfest venn-bar colours exactly.
INK, INK2, GRID = "#161616", "#4a4a46", "#e7e5dd"
BLUE, ORANGE, GREEN, RED = "#1f6fb2", "#e8720c", "#2f9e44", "#d1352b"
GOLD, PURPLE, TEAL = "#e0a82e", "#7b4ea3", "#1b9e9e"
COMMON, GAINED, LOST = "#115795", "#007D3E", "#E17224"  # native Oktoberfest gains/losses colours
TARGET_C, DECOY_C = "#2f9e44", "#8f8e88"  # target = green, decoy = neutral grey
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.edgecolor": "#c9c9c4",
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 110,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.labelcolor": INK2,
    }
)
FDR = 0.01


def fig_to_uri(fig, dpi=150):
    """Render a matplotlib figure to an inline PNG data URI and close it."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def svg_to_uri(path):
    """Read an SVG file and return it as an inline data URI, or "" if unreadable."""
    return "data:image/svg+xml;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def grid(ax, axis="both"):
    """Apply the report's house grid style to an axis."""
    ax.grid(axis=axis, color=GRID, lw=0.7)
    ax.set_axisbelow(True)


# ============================== scale guards ==============================
# This report is built inside a rescoring run, so every stage has to stay bounded when that run is a
# BULK experiment rather than a single-cell one: hundreds of raw files and tens of millions of PSMs
# instead of a few hundred thousand. The constants below cap what is parsed into memory, fitted,
# drawn and embedded. None of them changes what a panel shows on a small run.
CHUNK_ROWS = 200_000  # rows parsed per chunk -> parse memory is O(chunk), not O(file)
MAX_KDE_POINTS = 200_000  # violins fit a gaussian KDE, which costs O(n) per evaluation point
MAX_SCATTER_POINTS = 55_000  # points actually drawn in a density scatter (file size, not accuracy)
MAX_SPECTRA_FILES = 8  # mzML files opened for the spectra viewer
MAX_SVG_MB = 8.0  # a single native SVG above this is a scatter of millions of vector points

# Columns of rescore.tab, by how they are stored. Everything else becomes float32: on 10^7 PSMs the
# difference between float32 and the float64 default is gigabytes, and no panel resolves 7 digits.
F64_COLS = {"ScanNr", "Mass", "ExpMass"}  # exactness matters: scan numbers and precursor mass
CAT_COLS = {"filename", "Peptide"}  # few distinct values per run -> integer codes + levels
HASH_COLS = {"SpecId"}  # only ever joined on -> 64-bit hashes, never kept as text


class ReportTooLargeError(RuntimeError):
    """Raised when a run exceeds a guard rail, so the report is skipped deliberately, not attempted."""


# ============================== data loading ==============================
FDR_METHODS = ("percolator", "mokapot")
# Percolator and mokapot write the same set of files with different column names; the report reads both.
ID_COLS = ("PSMId", "SpecId")
SCORE_COLS = ("score", "mokapot score")
Q_COLS = ("q-value", "mokapot q-value")
PEPTIDE_COLS = ("peptide", "Peptide")


def _pick(header, names, path):
    """Index of the first of `names` present in `header`, so percolator and mokapot files both parse."""
    for name in names:
        if name in header:
            return header.index(name)
    raise KeyError(f"none of the columns {names} found in {path}")


def _hash_ids(values):
    """64-bit hashes of a sequence of identifier strings.

    PSM and peptide identifiers are only ever compared for equality (joining rescore.tab against the
    percolator/mokapot output, intersecting accepted sets) and never displayed, so they are stored as
    8-byte hashes instead of ~100-byte Python strings: on a bulk run that is the difference between a
    gigabyte-scale dict and a numpy array, and it turns every join into a vectorised operation.
    A collision (~1e-6 likely at 10^7 ids) would misplace one PSM in one diagnostic panel.
    """
    return np.fromiter((hash(v) for v in values), dtype=np.int64, count=len(values)).view(np.uint64)


def _lookup(keys, table_keys, table_values):
    """Vectorised dict-like lookup: table_values for each key, NaN where the key is absent."""
    out = np.full(keys.size, np.nan, dtype=np.float32)
    if keys.size == 0 or table_keys.size == 0:
        return out
    order = np.argsort(table_keys, kind="stable")
    sorted_keys = table_keys[order]
    pos = np.clip(np.searchsorted(sorted_keys, keys), 0, sorted_keys.size - 1)
    hit = sorted_keys[pos] == keys
    out[hit] = table_values[order][pos[hit]]
    return out


def _to_float(values, dtype):
    """Parse a chunk of string cells to floats, tolerating the odd empty or non-numeric one."""
    try:
        return np.array(values, dtype=dtype)
    except ValueError:  # one bad cell must cost that cell, not the whole column

        def one(v):
            try:
                return float(v)
            except ValueError:
                return np.nan

        return np.fromiter((one(v) for v in values), dtype=dtype, count=len(values))


def _column_chunk(col, values, code_map):
    if col in HASH_COLS:
        return _hash_ids(values)
    if col in CAT_COLS:
        return np.fromiter((code_map.setdefault(v, len(code_map)) for v in values), dtype=np.int32, count=len(values))
    return _to_float(values, np.float64 if col in F64_COLS else np.float32)


def find_paths(base):
    """Locate the FDR-method output dir, the spectra dir and the config of a run, from any of its dirs."""
    base = Path(base).resolve()
    cand = [base]
    cand += [base / "results" / m for m in FDR_METHODS]
    cand += [base / m for m in FDR_METHODS]
    hit = next(((c, m) for c in cand for m in FDR_METHODS if (c / f"rescore.{m}.psms.txt").exists()), None)
    if hit is None:  # last resort: walk the tree, which on a bulk run's output dir is not cheap
        for m in FDR_METHODS:
            found = next(base.rglob(f"rescore.{m}.psms.txt"), None)
            if found is not None:
                hit = (found.parent, m)
                break
    if hit is None:
        raise FileNotFoundError(f"could not find rescore.({'|'.join(FDR_METHODS)}).psms.txt under {base}")
    fdr_dir, method = hit
    # spectra dir + config: look near an `out` root
    out_root = base if (base / "spectra").exists() else fdr_dir.parent.parent
    spectra = next((d for d in [out_root / "spectra", base / "spectra"] if d.exists()), None)
    run_dir = out_root.parent
    config = next(
        (c for c in [run_dir / "config.json", out_root / "config.json", base / "config.json"] if c.exists()), None
    )
    results = fdr_dir.parent  # out/results (holds the CE-violin SVGs) is the FDR dir's parent
    return {
        "fdr_dir": fdr_dir,
        "method": method,
        "spectra": spectra,
        "config": config,
        "run_dir": run_dir,
        "out_root": out_root,
        "results_dir": out_root / "results" if (out_root / "results").exists() else results,
    }


def _chunks_pandas(path, header, cols, idx):
    """Chunks of raw cell values via pandas' C tokenizer.

    It parses numbers straight into arrays instead of building a Python object per cell, which on a
    bulk-size rescore.tab is several times faster and a fraction of the memory. Positional names plus
    index_col=False are what let it read a file whose trailing Proteins column carries extra tabs:
    the wanted columns are addressed by position, and everything past them is ignored rather than
    shifted into an index.
    """
    import pandas as pd

    as_text = {
        i for c, i in zip(cols, idx, strict=False) if c in CAT_COLS | HASH_COLS
    }  # ids stay strings, never numbers
    reader = pd.read_csv(
        path,
        sep="\t",
        names=range(len(header)),
        skiprows=1,
        usecols=idx,
        index_col=False,
        engine="c",
        chunksize=CHUNK_ROWS,
        dtype=dict.fromkeys(as_text, str),
        float_precision="round_trip",
    )  # exact: this path must match the fallback bit for bit
    for chunk in reader:
        yield {c: chunk[i] for c, i in zip(cols, idx, strict=False)}


def _chunks_python(path, cols, idx):
    """Chunks of raw cell values, parsed here. The fallback for a file pandas cannot tokenize."""
    mx = max(idx)
    buf = {c: [] for c in cols}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= mx:  # truncated line (e.g. an interrupted write): skip the row, keep the file
                continue
            for c, i in zip(cols, idx, strict=False):
                buf[c].append(p[i])
            if len(buf[cols[0]]) >= CHUNK_ROWS:
                yield buf
                buf = {c: [] for c in cols}
    if buf[cols[0]]:
        yield buf


def _read_columns(path, cols, idx, header, max_rows=0, log=None):
    """Assemble the columns at positions `idx` of a ragged tab file into arrays.

    Columns in CAT_COLS come back as int32 codes plus a `<col>_levels` array of the distinct strings,
    columns in HASH_COLS as 64-bit hashes, everything else as float32 (float64 for F64_COLS). The fast
    path is pandas' C tokenizer; if it cannot tokenize the file, the hand parser takes over and
    produces bit-identical arrays.

    :param path: the file to read
    :param cols: names to key the result on — they decide the dtype, see the module constants
    :param idx: column position in the file for each name
    :param header: the file's header line, already split
    :param max_rows: refuse files longer than this many rows (0 = no limit)
    :param log: optional progress sink, used to report a fallback to the slower parser
    :raises ReportTooLargeError: if the file holds more than max_rows rows
    :return: (dict of arrays, number of rows read)
    """

    def assemble(chunks):
        code_maps = {c: {} for c in cols if c in CAT_COLS}
        parts = {c: [] for c in cols}
        n = 0
        for raw in chunks:
            n += len(raw[cols[0]])
            if 0 < max_rows < n:
                raise ReportTooLargeError(
                    f"{Path(path).name} holds more than {max_rows:,} PSMs; raise the report option "
                    f"'max_psms' to build the report for a run this size"
                )
            for c in cols:
                parts[c].append(_column_chunk(c, raw[c], code_maps.get(c)))
        out = {}
        for c in cols:
            done = parts.pop(c)  # released column by column: holding chunks and result at once doubles the peak
            out[c] = np.concatenate(done) if done else _column_chunk(c, [], code_maps.get(c))
            del done
            if c in code_maps:
                out[f"{c}_levels"] = np.array(list(code_maps[c]), dtype=object)
        return out, n

    try:
        return assemble(_chunks_pandas(path, header, cols, idx))
    except ReportTooLargeError:
        raise
    except Exception as e:  # a file the C tokenizer chokes on is still readable by hand
        if log:
            log(f"[investigate]   fast parse of {Path(path).name} failed ({e}), falling back")
        return assemble(_chunks_python(path, cols, idx))


def read_tab(path, wanted, max_rows=0, log=None):
    """Read the wanted columns of rescore.tab that EXIST; see :py:func:`_read_columns` for the dtypes.

    Only the wanted columns are ever materialised: a bulk run's rescore.tab has tens of millions of
    rows and a hundred-odd feature columns, and reading it whole is not an option.

    :param path: the rescore.tab to read
    :param wanted: columns to keep, if present
    :param max_rows: refuse files longer than this many rows (0 = no limit)
    :param log: optional progress sink
    :raises ValueError: if none of the wanted columns is present
    :return: (dict of arrays, header, number of rows read)
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    cols = [c for c in wanted if c in header]
    if not cols:
        raise ValueError(f"none of the expected columns found in {path}")
    out, n = _read_columns(path, cols, [header.index(c) for c in cols], header, max_rows=max_rows, log=log)
    return out, header, n


def load_scores(path, key="id", log=None):
    """Read (id-hash, score, q-value) arrays from a percolator/mokapot psms or peptides output.

    `key` selects what the id hashes identify: the PSM ("id") or the peptide sequence ("peptide"),
    which is what the peptide-level files are keyed on. Rows whose score or q-value is unparseable
    come back as NaN, which places them outside every threshold the report applies.

    :param path: the psms/peptides output file to read
    :param key: "id" or "peptide"
    :param log: optional progress sink
    :return: (id hashes, scores, q-values)
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    idx = [
        _pick(header, PEPTIDE_COLS if key == "peptide" else ID_COLS, path),
        _pick(header, SCORE_COLS, path),
        _pick(header, Q_COLS, path),
    ]
    out, _ = _read_columns(path, ["SpecId", "score", "q"], idx, header, log=log)  # SpecId -> hashed id
    return out["SpecId"], out["score"], out["q"]


def ids_at_fdr(path, key="id", fdr=FDR):
    """Sorted unique id hashes accepted at `fdr` — set arithmetic on these is exact and vectorised."""
    ids, _, q = load_scores(path, key=key)
    return np.unique(ids[q <= fdr])


def parse_weights(path):
    """Percolator weights.csv -> (feature_names, mean normalized weight per feature).

    Layout after the 3 comment lines: (names, normalized, raw) repeated per CV bin.
    """
    lines = [ln.rstrip("\n") for ln in Path(path).read_text().splitlines() if not ln.startswith("#")]
    if not lines:
        return [], {}
    names = lines[0].split("\t")
    norm_rows = []
    for i in range(0, len(lines), 3):  # each triple: names, normalized, raw
        if i + 1 < len(lines):
            try:
                norm_rows.append([float(x) for x in lines[i + 1].split("\t")])
            except ValueError:
                pass
    if not norm_rows:
        return names, {}
    mean = np.mean([r for r in norm_rows if len(r) == len(names)], axis=0)
    return names, dict(zip(names, mean, strict=False))


def strip_pep(p):
    """Strip the flanking ``_.``/``._`` and the UNIMOD annotations from a percolator peptide."""
    return p[2:-2] if isinstance(p, str) and p.startswith("_.") and p.endswith("._") else p


def subsample(*arrays, k):
    """Deterministic uniform thinning to at most k rows, applied to every array in step.

    Used where a panel's message is a SHAPE (a density, a violin, a cloud) rather than a count: those
    are unchanged by thinning, while the cost of drawing or fitting them is not.
    """
    n = arrays[0].size
    if n <= k:
        return arrays if len(arrays) > 1 else arrays[0]
    sel = np.linspace(0, n - 1, k).astype(np.int64)
    out = tuple(a[sel] for a in arrays)
    return out if len(out) > 1 else out[0]


# ============================== load everything ==============================
def load_data(paths, max_psms=0, log=None):
    """Load rescore.tab and join the rescored / original scores and q-values onto every row.

    :param P: paths of the run, as returned by :py:func:`find_paths`
    :param max_psms: refuse runs with more PSMs than this (0 = no limit)
    :param log: optional progress sink
    :raises ValueError: if rescore.tab lacks the columns everything else is keyed on
    :return: dict of per-PSM arrays, plus `<col>_levels` for the coded columns and n / header
    """
    fdr_dir, method = paths["fdr_dir"], paths["method"]
    wanted = [
        "SpecId",
        "Label",
        "ScanNr",
        "filename",
        "ExpMass",
        "Mass",
        "RT",
        "iRT",
        "pred_RT",
        "abs_rt_diff",
        "collision_energy_aligned",
        "missedCleavages",
        "sequence_length",
        "KR",
        "delta_mass_ppm",
        "mean_ppm_error",
        "max_ppm_error",
        "log10_evalue",
        "pearson_corr",
        "spectral_angle",
        "spectral_angle_no_b1",
        "spectral_angle_noise_aware",
        "fraction_observed_and_predicted",
        "count_observed_and_predicted",
        "count_predicted",
        "Peptide",
        "Charge1",
        "Charge2",
        "Charge3",
        "Charge4",
        "Charge5",
        "Charge6",
    ]
    t, header, n = read_tab(fdr_dir / "rescore.tab", wanted, max_rows=max_psms, log=log)
    for required in ("SpecId", "Label"):
        if required not in t:
            raise ValueError(f"rescore.tab has no {required} column, cannot build the report")
    # a truncated or garbled row has no usable label; drop it rather than let it count as a third class
    usable = np.isin(t["Label"], (-1, 1))
    if not usable.all():
        levels = {k: v for k, v in t.items() if k.endswith("_levels")}  # per-column, not per-row
        t = {k: v[usable] for k, v in t.items() if k not in levels}
        t.update(levels)
        if log:
            log(f"[investigate]   dropped {int((~usable).sum()):,} row(s) without a valid Label")
        n = int(usable.sum())
    charge = np.zeros(n, np.int8)
    for c in range(1, 7):
        key = f"Charge{c}"
        if key in t:
            charge[t[key] == 1] = c
    label = t["Label"].astype(np.int8)

    def scores_of(prefix):
        """(score, q) per row from this run's `prefix` output; q stays target-only, as every panel means it."""
        psms = fdr_dir / f"{prefix}.{method}.psms.txt"
        if not psms.exists():
            return np.full(n, np.nan, np.float32), np.full(n, np.nan, np.float32)
        tid, tscore, tq = load_scores(psms)
        decoys = fdr_dir / f"{prefix}.{method}.decoy.psms.txt"
        did, dscore = load_scores(decoys)[:2] if decoys.exists() else (np.empty(0, np.uint64), np.empty(0, np.float32))
        # target and decoy ids are disjoint, so one table serves both label classes
        score = _lookup(t["SpecId"], np.concatenate([tid, did]), np.concatenate([tscore, dscore]))
        return score, _lookup(t["SpecId"], tid, tq)

    score, q = scores_of("rescore")
    oscore, oq = scores_of("original")
    t.update(
        {"charge": charge, "label": label, "score": score, "q": q, "oscore": oscore, "oq": oq, "header": header, "n": n}
    )
    return t


# ============================== FIGURES ==============================
RUN_CTX = {}  # {"suffix": "<engine> · <model> · <run>"} set by main() -> adaptive figure subtitles


def _caption(ax, text):
    ax.set_title(text, fontsize=10.5, loc="left", color=INK, fontweight="bold", pad=8)


def _sup(fig, main, y=0.998):
    fig.suptitle(main, fontsize=14, fontweight="bold", y=y)
    sub = RUN_CTX.get("suffix", "")
    if sub:
        gap = 0.42 / max(fig.get_figheight(), 3)  # ~constant visual gap; keeps short figs from colliding
        fig.text(0.5, y - gap, sub, ha="center", va="top", fontsize=10, color=INK2, style="italic")


def density_scatter(ax, x, y, title, xlabel, clip=(1, 99), ylab="spectral angle", ylim=(0, 1.02), trend=True):
    """Native-Oktoberfest-style scatter, coloured by local 2-D density.

    Warm cmap, log scale, dense points drawn on top, so sparse vs dense regions stay legible.
    """
    from scipy.stats import spearmanr

    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 20:
        ax.axis("off")
        return
    lo, hi = np.percentile(x, clip)
    k = (x >= lo) & (x <= hi)
    x, y = x[k], y[k]
    n_true = x.size
    nb = 64
    hist, xe, ye = np.histogram2d(x, y, bins=nb)
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, nb - 1)
    iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, nb - 1)
    dens = hist[ix, iy]
    # trend + correlation on the FULL data
    r, _ = spearmanr(x, y)
    mx, my = [], []
    if trend:
        bins = np.linspace(x.min(), x.max(), 12)
        bi = np.clip(np.digitize(x, bins), 1, len(bins) - 1)
        for b in range(1, len(bins)):
            s = bi == b
            if s.sum() > 20:
                mx.append(x[s].mean())
                my.append(np.median(y[s]))
    # scatter (subsample for file size, keeping the density ordering so dense stays on top)
    o = np.argsort(dens, kind="stable")
    xs, ys, ds = x[o], y[o], dens[o]
    xs, ys, ds = subsample(xs, ys, ds, k=MAX_SCATTER_POINTS)  # thins the drawing, not the statistics
    sc = ax.scatter(xs, ys, c=np.log1p(ds), cmap="YlOrRd", s=5, linewidths=0, rasterized=True)
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_ticks([])
    cb.ax.set_ylabel("point density", fontsize=8, color=INK2)
    cb.ax.text(0.5, 1.01, "dense", transform=cb.ax.transAxes, ha="center", va="bottom", fontsize=7, color=INK2)
    cb.ax.text(0.5, -0.01, "sparse", transform=cb.ax.transAxes, ha="center", va="top", fontsize=7, color=INK2)
    if trend:
        ax.plot(mx, my, "-o", color=INK, lw=1.8, ms=3, label="binned median")
        ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.text(
        0.97,
        0.06,
        f"Spearman r={r:+.2f}\nn={n_true:,}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=INK,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#cccccc", "alpha": 0.85},
    )
    _caption(ax, title)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlabel)
    grid(ax)


def _violins(ax, groups, positions, color, width=0.75):
    # the KDE behind a violin costs O(n) per evaluation point; its shape is settled long before 200k
    v = ax.violinplot(
        [subsample(g, k=MAX_KDE_POINTS) for g in groups], positions=positions, widths=width, showextrema=False
    )
    for b in v["bodies"]:
        b.set_facecolor(color)
        b.set_edgecolor(INK2)
        b.set_alpha(0.55)
        b.set_linewidth(1)
    for g, p in zip(groups, positions, strict=False):
        med = np.median(g)
        ax.hlines(med, p - width / 3, p + width / 3, color=INK, lw=1.8)
        ax.text(p, 1.03, f"{med:.2f}", ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")


def _split_violin(ax, left, right, col_left, col_right, pos=0.0, width=0.9):
    """One violin split in two: left half = `left` distribution, right half = `right`."""
    for data, side, col in [(left, "left", col_left), (right, "right", col_right)]:
        v = ax.violinplot([subsample(data, k=MAX_KDE_POINTS)], positions=[pos], widths=width, showextrema=False)
        b = v["bodies"][0]
        verts = b.get_paths()[0].vertices
        if side == "left":
            verts[:, 0] = np.clip(verts[:, 0], -np.inf, pos)
        else:
            verts[:, 0] = np.clip(verts[:, 0], pos, np.inf)
        b.set_facecolor(col)
        b.set_edgecolor(INK2)
        b.set_alpha(0.6)
        b.set_linewidth(1)
    ax.hlines(np.median(left), pos - width / 2, pos, color=INK, lw=2)
    ax.hlines(np.median(right), pos, pos + width / 2, color=INK, lw=2)
    ax.text(
        pos - width / 4,
        1.03,
        f"{np.median(left):.2f}",
        ha="center",
        va="bottom",
        fontsize=9,
        color=INK,
        fontweight="bold",
    )
    ax.text(
        pos + width / 4,
        1.03,
        f"{np.median(right):.2f}",
        ha="center",
        va="bottom",
        fontsize=9,
        color=INK,
        fontweight="bold",
    )


def fig_sa_features(data, paths):
    """Build the "spectral angle vs peptide and spectrum properties" panel."""
    if "spectral_angle" not in data:
        return None  # a feature-ablation run can drop it; every SA panel is then vacuous
    acc = (data["label"] == 1) & (data["q"] <= FDR)
    sa = data["spectral_angle"]
    saa = sa[acc]
    ch = data["charge"][acc]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
    panels = axes.ravel()
    # A0: SA by precursor charge — split violin (2+ | 3+) when both present, else side-by-side violins
    c2, c3 = saa[ch == 2], saa[ch == 3]
    if len(c2) > 30 and len(c3) > 30:
        _split_violin(panels[0], c2, c3, BLUE, ORANGE)
        panels[0].set_xticks([])
        panels[0].set_xlim(-0.7, 0.7)
        panels[0].legend(
            handles=[
                Patch(facecolor=BLUE, alpha=0.6, edgecolor=INK2, label=f"charge 2+  (n={len(c2):,})"),
                Patch(facecolor=ORANGE, alpha=0.6, edgecolor=INK2, label=f"charge 3+  (n={len(c3):,})"),
            ],
            frameon=False,
            fontsize=9,
            loc="lower center",
        )
        _caption(panels[0], "SA by precursor charge  (split violin: 2+ | 3+)")
    else:
        chs = [c for c in [1, 2, 3, 4, 5, 6] if (ch == c).sum() > 30]
        _violins(panels[0], [saa[ch == c] for c in chs], list(range(len(chs))), BLUE)
        panels[0].set_xticks(range(len(chs)))
        panels[0].set_xticklabels([f"{c}+\n(n={int((ch == c).sum()):,})" for c in chs])
        panels[0].set_xlabel("precursor charge")
        _caption(panels[0], "SA by precursor charge  (violins)")
    panels[0].set_ylim(0, 1.12)
    panels[0].set_ylabel("spectral angle")
    grid(panels[0], "y")
    # A1: SA by missed cleavages — violins
    mc = np.clip(data["missedCleavages"][acc], 0, 2).astype(int)
    mcg = [saa[mc == k] for k in [0, 1, 2]]
    _violins(panels[1], mcg, [0, 1, 2], GOLD)
    panels[1].set_xticks([0, 1, 2])
    panels[1].set_xticklabels([f"{'2+' if k == 2 else k}\n(n={len(mcg[k]):,})" for k in [0, 1, 2]])
    panels[1].set_ylim(0, 1.12)
    panels[1].set_ylabel("spectral angle")
    panels[1].set_xlabel("missed cleavages")
    _caption(panels[1], "SA by missed cleavages  (violins)")
    grid(panels[1], "y")
    # A2..A3: the two continuous covariates that actually move SA. Deliberately NOT plotted, because each
    # was flat or duplicated a panel we keep: retention time (Spearman ~+0.05 — no effect), matched-fragment
    # COUNT and peptide LENGTH (collinear with fragment coverage / precursor mass but weaker), and SA-vs-RT-
    # residual (weak cross-term; the RT residual gets its own target-vs-decoy panel under Calibration).
    scatter_panels = [
        ("fraction_observed_and_predicted", "fragment coverage", "fraction obs & pred"),
        ("Mass", "precursor mass", "precursor mass (Da)"),
    ]
    ai = 2
    for col, title, lab in scatter_panels:
        if ai >= 4:
            break
        if col not in data or not np.isfinite(data[col][acc]).any():
            continue
        density_scatter(panels[ai], data[col][acc], saa, title, lab)
        ai += 1
    for j in range(ai, 4):
        panels[j].axis("off")
    _sup(fig, "SA diagnostics — spectral angle vs features (confident targets, q≤1%)")
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=3.2, w_pad=2.2)
    cap = (
        "Confident target PSMs (q≤1%), showing the four properties that measurably move the spectral angle. "
        "Violin panels: shaded = SA distribution, black bar = median (value above); SA-by-charge is a split "
        "violin — left half 2+, right half 3+. Scatter panels plot every PSM coloured by local point density "
        "(warm colourbar = dense), with a black binned-median line; Spearman r and n annotated, x clipped to "
        "[p1,p99]. Reading: higher charge, missed cleavages and heavier precursors all lower SA, while fragment "
        "coverage raises it — i.e. SA degrades exactly where the fragmentation model is least constrained."
    )
    return "sa", "SA diagnostics", cap, fig_to_uri(fig)


def fig_headroom(data, paths):
    """Build the "rescoring headroom" panel: rejected targets against decoys."""
    if "spectral_angle" not in data:
        return None
    sa, label, q = data["spectral_angle"], data["label"], data["q"]
    acc = (label == 1) & (q <= FDR)
    rej = (label == 1) & (q > FDR)
    isd = label == -1
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8))
    bins = np.linspace(0, 1, 60)
    # (a)
    ax = axes[0]
    for m, col, lbl in [
        (acc, GREEN, f"accepted targets (q≤1%), n={int(acc.sum()):,}"),
        (rej, ORANGE, f"rejected targets (q>1%), n={int(rej.sum()):,}"),
        (isd, INK2, f"decoys, n={int(isd.sum()):,}"),
    ]:
        ax.hist(sa[m], bins=bins, density=True, histtype="step", lw=2, color=col, label=lbl)
    ax.set_xlabel("spectral angle")
    ax.set_ylabel("density (area=1)")
    ax.legend(frameon=False, fontsize=9)
    _caption(ax, "(a) SA: accepted vs rejected targets vs decoys")
    grid(ax)
    # (b) — the quantitative version of (a): would an SA cut recover anything the decoys don't also pass?
    ax = axes[1]
    thr = [0.5, 0.6, 0.7, 0.8]
    rc = [int((rej & (sa > th)).sum()) for th in thr]
    dc = [int((isd & (sa > th)).sum()) for th in thr]
    x = np.arange(len(thr))
    w = 0.38
    ax.bar(x - w / 2, rc, w, color=ORANGE, label="rejected targets (q>1%)")
    ax.bar(x + w / 2, dc, w, color=INK2, label="decoys")
    top = max(max(rc), max(dc), 1)
    for i, (r, d) in enumerate(zip(rc, dc, strict=False)):
        ax.text(
            i,
            max(r, d) + top * 0.02,
            f"net {r - d:+,}",
            ha="center",
            va="bottom",
            fontsize=10,
            color=GREEN if r - d > 0 else RED,
            fontweight="bold",
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"SA>{th}" for th in thr])
    ax.set_ylabel("PSM count")
    ax.set_ylim(0, max(max(rc), max(dc), 1) * 1.2)
    ax.legend(frameon=False, fontsize=9)
    _caption(ax, "(b) high-SA rejected targets vs high-SA decoys")
    grid(ax, "y")
    _sup(fig, "Rescoring headroom — target PSMs NOT accepted at 1% FDR", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    # The caption states how to READ the net numbers, not what they come out as: the sign is run-dependent
    # (it is positive on some search-engine/model combinations and negative on others), so asserting a
    # direction here would make the report wrong on half the runs it is generated for.
    cap = (
        "Do we leave good identifications behind? (a) the SA distribution of REJECTED targets, compared against "
        "the DECOYS: where the two overlap, the rejects are indistinguishable from noise; any excess of rejected "
        "targets over decoys at high SA is the part that could still be real. (b) the same comparison in counts. "
        "Decoys estimate how many false PSMs an SA cut would also let through, so net = (rejected − decoy) "
        "approximates the true targets such a cut could recover: a large positive net means real headroom that "
        "a better spectral model could reach, while a net near zero or negative means an SA cut would buy nothing "
        "the FDR has not already priced in."
    )
    return "headroom", "Rescoring headroom", cap, fig_to_uri(fig)


def fig_weights(data, paths):
    """Build the percolator/mokapot feature-weight panel, with vs without Prosit."""
    # weights.csv is a percolator artefact; mokapot writes none, so this panel is percolator-only
    resc = paths["fdr_dir"] / "rescore.percolator.weights.csv"
    orig = paths["fdr_dir"] / "original.percolator.weights.csv"
    if not resc.exists():
        return None
    rn, rw = parse_weights(resc)
    on, ow = parse_weights(orig) if orig.exists() else ([], {})
    two = bool(ow)
    fig, axes = plt.subplots(1, 2 if two else 1, figsize=(14.5 if two else 8, 9.2), squeeze=False)

    def panel(ax, w, title, added=None):
        if added is None:
            added = set()
        items = [(k, v) for k, v in w.items() if abs(v) > 1e-4]  # drop exact-zero (unused) features
        items = sorted(items, key=lambda kv: abs(kv[1]), reverse=True)[:22][::-1]
        names = [k for k, _ in items]
        vals = [v for _, v in items]
        cols = [
            (GREEN if names[i] in added else BLUE) if vals[i] >= 0 else (RED if names[i] not in added else GOLD)
            for i in range(len(names))
        ]
        ax.barh(range(len(names)), vals, color=cols)
        ax.axvline(0, color=INK, lw=1)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=8.5)
        ax.set_xlabel("normalized SVM weight  (mean over CV bins)")
        _caption(ax, title)
        grid(ax, "x")

    added = set(rn) - set(on)
    panel(axes[0][0], rw, "+Prosit (all features)", added)
    if two:
        panel(axes[0][1], ow, "no Prosit (search features only)")
    hs = [
        Patch(color=BLUE, label="search feature, +weight (→ target)"),
        Patch(color=RED, label="search feature, −weight"),
        Patch(color=GREEN, label="Prosit-added feature, +weight"),
        Patch(color=GOLD, label="Prosit-added feature, −weight"),
    ]
    fig.legend(handles=hs, loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    _sup(fig, "Percolator feature weights — what drives the score", y=0.995)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    cap = (
        "Mean normalized Percolator SVM weight per feature (averaged over cross-validation bins), top 22 by "
        "|weight| (exact-zero/unused features dropped). Positive (right) pushes a PSM toward TARGET. Left = the "
        "+Prosit model (green = features Prosit adds: spectral_angle*, pearson_corr, RT, fragment counts…); right = "
        "the no-Prosit model. NB Oktoberfest names the primary search score 'andromeda'/'m0' even for MSFragger."
    )
    return "weights", "Percolator feature weights", cap, fig_to_uri(fig)


def _qvals(path):
    """Sorted q-values of the target rows in a percolator/mokapot output file."""
    return np.sort(load_scores(path)[2])


def fig_yield(data, paths):
    """Identifications accepted as a function of the FDR threshold, +Prosit vs no-Prosit.

    The standard way this comparison is reported in the rescoring literature: one curve per model, so the
    gain is visible at EVERY threshold rather than only at the single 1% operating point.
    """
    lv_have = []
    for lv, name in [("psms", "PSMs"), ("peptides", "Peptides")]:
        r = paths["fdr_dir"] / f"rescore.{paths['method']}.{lv}.txt"
        o = paths["fdr_dir"] / f"original.{paths['method']}.{lv}.txt"
        if r.exists():
            lv_have.append((name, _qvals(r), _qvals(o) if o.exists() else None))
    if not lv_have:
        return None
    ts = np.linspace(0, 0.05, 501)
    fig, axes = plt.subplots(1, len(lv_have), figsize=(7.1 * len(lv_have), 5.6), squeeze=False)
    for i, (name, qr, qo) in enumerate(lv_have):
        ax = axes[0][i]
        yr = np.searchsorted(qr, ts, side="right")
        if qo is not None:
            yo = np.searchsorted(qo, ts, side="right")
            ax.fill_between(ts, yo, yr, color=GAINED, alpha=0.16, lw=0, label="gain from Prosit")
            ax.plot(ts, yo, color=COMMON, lw=2.2, label="no Prosit (search features only)")
        ax.plot(ts, yr, color=GAINED, lw=2.2, label="+Prosit (all features)")
        ax.axvline(FDR, color=RED, ls="--", lw=1.4, label="1% FDR")
        nr = int(np.searchsorted(qr, FDR, side="right"))
        ax.plot([FDR], [nr], "o", color=GAINED, ms=6, zorder=5)
        txt = f"{nr:,}"
        if qo is not None:
            no = int(np.searchsorted(qo, FDR, side="right"))
            ax.plot([FDR], [no], "o", color=COMMON, ms=6, zorder=5)
            pct = 100 * (nr - no) / no if no else float("nan")
            txt = f"{nr:,}  vs  {no:,}\n{nr - no:+,}  ({pct:+.1f}%)"
        ax.annotate(
            txt,
            (FDR, nr),
            textcoords="offset points",
            xytext=(12, -6),
            fontsize=10,
            color=INK,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.32", "fc": "white", "ec": "#cccccc", "alpha": 0.9},
        )
        ax.set_xlim(0, 0.05)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("FDR threshold (q-value)")
        ax.set_ylabel(f"accepted target {name} (cumulative)")
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        _caption(ax, f"{name} accepted vs FDR threshold")
        grid(ax)
    _sup(fig, "Identification yield vs FDR — with and without Prosit", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    fdr = paths["method"].capitalize()
    cap = (
        f"Number of target identifications accepted as the FDR threshold is relaxed, for {fdr} WITH Prosit "
        "predictions (green) and WITHOUT them (blue); the shaded band between the curves is the gain. The dashed "
        "line marks the 1% operating point, annotated with both counts, the absolute and the relative gain. "
        "Reading: a curve that sits above the other ACROSS the whole range is a genuine discrimination "
        "improvement, not an artefact of where the threshold happens to be drawn."
    )
    return "yield", "Identification yield vs FDR", cap, fig_to_uri(fig)


def fig_movement(data, paths):
    """Build the per-PSM score-movement panel, coloured by what was gained, lost or kept."""
    label, score, oscore = data["label"], data["score"], data["oscore"]
    m = (label == 1) & np.isfinite(score) & np.isfinite(oscore)
    if m.sum() < 50:
        return None
    # acceptance at 1% FDR is what the per-row q-values already say: no id set membership needed
    # (that was a Python loop over every PSM, plus two sets of every accepted id, at bulk scale)
    in_r = data["q"] <= FDR
    in_o = data["oq"] <= FDR
    gained = m & in_r & ~in_o
    lost = m & ~in_r & in_o
    common = m & in_r & in_o
    neither = m & ~in_r & ~in_o
    fig, ax = plt.subplots(figsize=(8.6, 8))
    rng = np.random.default_rng(0)

    def sub(mask, k=8000):
        idx = np.where(mask)[0]
        return rng.choice(idx, size=min(k, idx.size), replace=False) if idx.size else idx

    grey = sub(neither, 6000)  # one draw: x and y must come from the same PSMs
    ax.scatter(
        oscore[grey],
        score[grey],
        s=4,
        color="#d0d0cc",
        alpha=0.4,
        linewidths=0,
        label=f"neither, n={int(neither.sum()):,}",
    )
    for mask, col, lbl in [
        (common, COMMON, f"kept (both), n={int(common.sum()):,}"),
        (gained, GAINED, f"gained by Prosit, n={int(gained.sum()):,}"),
        (lost, LOST, f"lost by Prosit, n={int(lost.sum()):,}"),
    ]:
        s = sub(mask)
        ax.scatter(oscore[s], score[s], s=6, color=col, alpha=0.5, linewidths=0, label=lbl)
    lim = [np.nanpercentile(oscore[m], 0.5), np.nanpercentile(oscore[m], 99.5)]
    ax.plot(lim, lim, color=INK, lw=0.9, ls="--", label="no change")
    ax.set_xlabel("original score (no Prosit)")
    ax.set_ylabel("rescored score (+Prosit)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    _caption(ax, "Per-PSM score movement: original → +Prosit")
    grid(ax)
    _sup(fig, "What Prosit moved — original vs rescored score", y=0.99)
    cap = (
        f"Each point is a target PSM: x = its {paths['method'].capitalize()} score without Prosit, y = with Prosit. Colour = acceptance "
        "at 1% FDR — green: accepted only after Prosit (gained), orange: accepted only before (lost), dark-blue: "
        "accepted both, grey: neither. Points above the dashed y=x line were pushed up by Prosit."
    )
    return "movement", "Rescore-vs-original movement", cap, fig_to_uri(fig)


def _distinct_per_file(file_codes, pep_codes, mask, n_files):
    """Distinct peptides per raw file among the masked rows, without materialising any per-file mask."""
    pairs = np.stack([file_codes[mask].astype(np.int64), pep_codes[mask].astype(np.int64)], axis=1)
    if pairs.size == 0:
        return np.zeros(n_files, dtype=np.int64)
    return np.bincount(np.unique(pairs, axis=0)[:, 0], minlength=n_files)


def fig_per_file(data, paths):
    """Identifications PER RAW FILE, +Prosit vs no-Prosit.

    Every other panel in this report pools all raw files together, which hides the question that matters for
    single-cell / low-input runs: does rescoring lift every file, or only rescue a few good ones?
    """
    if "filename" not in data or not np.isfinite(data["oq"]).any():
        return None
    codes, n_files = data["filename"], data["filename_levels"].size
    if n_files < 4:
        return None  # a per-file view needs enough files to be a distribution rather than a list
    acc_r = (data["label"] == 1) & (data["q"] <= FDR)
    acc_o = (data["label"] == 1) & (data["oq"] <= FDR)
    # counted on the integer peptide codes: identical to counting the strings, but a run with hundreds
    # of files no longer costs one full-length mask per file (that loop is quadratic in the run size)
    if "Peptide" in data:
        unit = "distinct peptides"
        nr, no = (
            _distinct_per_file(codes, data["Peptide"], acc_r, n_files),
            _distinct_per_file(codes, data["Peptide"], acc_o, n_files),
        )
    else:
        unit = "PSMs"
        nr, no = np.bincount(codes[acc_r], minlength=n_files), np.bincount(codes[acc_o], minlength=n_files)
    nr = nr.astype(float)
    no = no.astype(float)
    order = np.argsort(-nr)
    x = np.arange(n_files)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    # (a) the per-file yield curve — dynamic range across files, and the lift on top of it
    ax = axes[0]
    ax.fill_between(x, no[order], nr[order], color=GAINED, alpha=0.18, lw=0, label="gain from Prosit")
    ax.plot(x, no[order], color=COMMON, lw=1.8, label="no Prosit")
    ax.plot(x, nr[order], color=GAINED, lw=1.8, label="+Prosit")
    ax.set_xlabel(f"raw file, ranked by +Prosit yield  (n={n_files} files)")
    ax.set_ylabel(f"{unit} @ 1% FDR")
    ax.set_xlim(0, n_files - 1)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9)
    _caption(ax, f"(a) {unit} per raw file  ·  median {np.median(nr):,.0f} vs {np.median(no):,.0f}")
    grid(ax)
    # (b) does the lift depend on how well the file ran?
    ax = axes[1]
    ok = no > 0
    pct = 100 * (nr[ok] - no[ok]) / no[ok]
    ax.scatter(no[ok], pct, s=26, color=GAINED, alpha=0.75, linewidths=0)
    ax.axhline(0, color=INK, lw=1)
    med = float(np.median(pct))
    ax.axhline(med, color=ORANGE, ls="--", lw=1.5, label=f"median {med:+.1f}%")
    n_up = int((nr[ok] > no[ok]).sum())
    n_dn = int((nr[ok] < no[ok]).sum())
    from scipy.stats import spearmanr

    # a rank correlation needs both axes to vary; with every file at the same yield or the same gain
    # there is no trend to report, and scipy would return NaN
    varies = np.unique(no[ok]).size > 1 and np.unique(pct).size > 1
    corr = f"Spearman r={spearmanr(no[ok], pct)[0]:+.2f}" if varies else "Spearman r: n/a (no spread)"
    ax.text(
        0.97,
        0.94,
        f"gained in {n_up}/{int(ok.sum())} files · lost in {n_dn}\n{corr}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color=INK,
        bbox={"boxstyle": "round,pad=0.32", "fc": "white", "ec": "#cccccc", "alpha": 0.9},
    )
    ax.set_xlabel(f"{unit} @ 1% FDR without Prosit  (per raw file)")
    ax.set_ylabel("relative gain from Prosit (%)")
    lo = min(0.0, float(pct.min()))
    hi = max(float(pct.max()), lo + 1e-6)  # never a zero-height axis
    ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.18 * (hi - lo))
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    _caption(ax, "(b) relative gain vs how well the file ran")
    grid(ax)
    _sup(fig, "Per-raw-file identifications — is the gain uniform?", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    cap = (
        f"Every other panel pools all raw files; this one splits them. (a) {unit} accepted at 1% FDR per raw file, "
        "files ranked by their +Prosit yield, with vs without Prosit — the spread along x is the run-to-run "
        "variability of the experiment, the shaded band is what rescoring adds on top of it. (b) the same gain "
        "expressed relative to each file's own no-Prosit yield, plotted against that yield: a flat cloud means "
        "rescoring helps every file by roughly the same factor, a downward trend (negative Spearman r) means "
        "weak files benefit disproportionately. NB counts are distinct peptide sequences among PSMs accepted at "
        f"1% PSM-level FDR, which is a per-file proxy — {paths['method'].capitalize()}'s own peptide-level FDR is estimated globally, "
        "not per file, so these do not sum to the peptide count in the summary."
    )
    return "perfile", "Per-raw-file identifications", cap, fig_to_uri(fig)


def fig_calibration(data, paths):
    """Build the retention-time and mass-accuracy calibration panel."""
    label, q = data["label"], data["q"]
    acc = (label == 1) & (q <= FDR)
    dec = label == -1
    have_ard = "abs_rt_diff" in data
    ppm_cols = [
        (c, lab)
        for c, lab in [
            ("delta_mass_ppm", "precursor mass error (ppm)"),
            ("mean_ppm_error", "fragment mean ppm error (ppm)"),
        ]
        if c in data
    ]
    # NB no pooled observed-vs-predicted iRT scatter here: RT alignment is fitted PER RAW FILE, so pooling
    # all raws bends the band into a curve that says more about the pooling than about the model. The honest
    # per-raw straight-line fits are the iRT-vs-RT SVGs in the per-raw gallery; the pooling-safe summary is
    # the abs_rt_diff residual below.
    panels = []  # list of draw-callables
    if have_ard:

        def _ard(ax):
            va = data["abs_rt_diff"][acc]
            vd = data["abs_rt_diff"][dec]
            va = va[np.isfinite(va)]
            vd = vd[np.isfinite(vd)]
            hi = np.percentile(np.concatenate([va, vd]), 99)
            b = np.linspace(0, hi, 60)
            ax.hist(
                va,
                bins=b,
                density=True,
                histtype="step",
                lw=2,
                color=TARGET_C,
                label=f"accepted targets, n={va.size:,}",
            )
            ax.hist(vd, bins=b, density=True, histtype="step", lw=2, color=DECOY_C, label=f"decoys, n={vd.size:,}")
            ax.axvline(np.median(va), color=TARGET_C, ls=":", lw=1.3)
            ax.axvline(np.median(vd), color=DECOY_C, ls=":", lw=1.3)
            ax.set_xlabel("|aligned RT − predicted RT|  (abs_rt_diff)")
            ax.set_ylabel("density (area=1)")
            ax.legend(frameon=False, fontsize=8.5)
            _caption(ax, f"RT residual (pooling-safe)  median tgt={np.median(va):.2f} vs decoy={np.median(vd):.2f}")
            grid(ax)

        panels.append(_ard)
    for c, lab in ppm_cols:

        def _ppm(ax, c=c, lab=lab):
            va = data[c][acc]
            vd = data[c][dec]
            va = va[np.isfinite(va)]
            vd = vd[np.isfinite(vd)]
            # robust MAD window around the target median: zooms to the monoisotopic core, so isotope-error
            # satellites (±1 Da ≈ hundreds of ppm) and broad decoys spill beyond instead of crushing it.
            med = np.median(va)
            mad = np.median(np.abs(va - med))
            half = 6 * mad if mad > 1e-6 else max(np.percentile(va, 97) - med, 1.0)
            lo, hi = max(va.min(), med - half), med + half
            if hi - lo < 1e-6:
                hi = lo + 1
            b = np.linspace(lo, hi, 60)
            ax.hist(
                va,
                bins=b,
                range=(lo, hi),
                density=True,
                histtype="step",
                lw=2,
                color=TARGET_C,
                label=f"accepted targets, n={va.size:,}",
            )
            ax.hist(
                vd,
                bins=b,
                range=(lo, hi),
                density=True,
                histtype="step",
                lw=2,
                color=DECOY_C,
                label=f"decoys, n={vd.size:,}",
            )
            ax.set_xlim(lo, hi)
            ax.set_xlabel(lab + "  [target core]")
            ax.set_ylabel("density (area=1)")
            ax.legend(frameon=False, fontsize=8.5)
            _caption(ax, lab.split(" (")[0])
            grid(ax)

        panels.append(_ppm)
    if not panels:
        return None
    nrow, ncol = (2, 2) if len(panels) > 2 else (1, len(panels))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 5.2 * nrow), squeeze=False)
    flat = axes.ravel()
    for i, draw in enumerate(panels):
        draw(flat[i])
    for j in range(len(panels), len(flat)):
        flat[j].axis("off")
    _sup(fig, "Calibration — retention time & mass accuracy", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94 if nrow > 1 else 0.9), h_pad=2.5, w_pad=2.5)
    cap = (
        "Each panel contrasts accepted targets (green) with decoys (grey) — the separation between the two IS the "
        f"discriminative power that feature contributes to {paths['method'].capitalize()}. RT panel: the pooling-safe per-PSM aligned RT "
        "residual |aligned RT − predicted RT| (targets tight, decoys broad; medians dotted). The per-raw "
        "observed-vs-predicted iRT fits are in the gallery below, where the per-file alignment is not pooled away. "
        "ppm panels: mass error, zoomed to the TARGET core via a robust MAD window (isotope-error PSMs spill beyond)."
    )
    return "calibration", "RT & mass-error calibration", cap, fig_to_uri(fig)


# ============================== run parameters ==============================
# Everything the report treats as a property of the EXPERIMENT rather than of the code — the isobaric
# label, the fragment match tolerance, the fragmentation method, the prediction server — is read from
# the run's own config.json here and passed around as one dict. Nothing downstream may hard-code such
# a value: the same report has to be right for a label-free bulk run and for a TMT single-cell one,
# and a constant that happens to fit the run in front of us silently misdescribes every other run.
# Where the config is silent, the value falls back to the same default OKTOBERFEST would have used
# (see :py:class:`oktoberfest.utils.config.Config`) and is MARKED as a default in the parameter
# table, so a default is never displayed as though the user had chosen it.

#: Oktoberfest's own defaults for the config keys the parameter table shows, by dotted path.
CFG_DEFAULTS = {
    "prediction_server": "koina.wilhelmlab.org:443",
    "ssl": True,
    "tag": "",
    "fragmentation_method": "HCD",
    "ion_types": "yb",
    "p_window": 0.0,
    "matching_method": "nearest",
    "matching_method_params": {},
    "fdr_estimation_method": "percolator",
    "allFeatures": False,
    "add_feature_cols": "none",
    "regressionMethod": "spline",
    "quantification": False,
    "inputs.search_results_type": "maxquant",
    "inputs.spectra_type": "raw",
    "ce_alignment_options.ce_range": [18, 50],
    "ce_alignment_options.use_ransac_model": False,
    "fastaDigestOptions.digestion": "full",
    "fastaDigestOptions.enzyme": "trypsin",
    "fastaDigestOptions.missedCleavages": 2,
    "fastaDigestOptions.minLength": 7,
    "fastaDigestOptions.maxLength": 60,
    "fastaDigestOptions.db": "concat",
}

#: Isobaric reporter-ion regions, keyed by the ``tag`` values Oktoberfest accepts (the keys of
#: spectrum_fundamentals' ``TMT_MODS``). Reporter ions are fragments of the LABEL, not of the peptide,
#: so the spectra viewer hides them — but only for a run that actually carried that label: in a
#: label-free (bulk) run the same region holds genuine immonium and low-mass fragment ions, and
#: blanking it would misrepresent the spectrum. Each window brackets that plex's reporter masses with
#: ~0.6 Da of margin: TMT 126.128–131.145, TMTpro 126.128–134.155, iTRAQ4 114.111–117.115,
#: iTRAQ8 113.108–121.122.
REPORTER_WINDOWS = {"tmt": (125.5, 132.0), "tmtpro": (125.5, 135.0), "itraq4": (113.5, 118.0), "itraq8": (112.0, 122.0)}
LABEL_NAMES = {
    "tmt": "TMT (2/6/10/11-plex)",
    "tmtpro": "TMTpro (16/18-plex)",
    "itraq4": "iTRAQ 4-plex",
    "itraq8": "iTRAQ 8-plex",
}

#: Fallback fragment tolerances per mass analyzer, mirroring
#: :py:func:`spectrum_fundamentals.fragments.get_min_max_mass` — i.e. what the RUN itself annotated
#: peaks with whenever the config set no explicit massTolerance/unitMassTolerance pair.
ANALYZER_TOL = {"FTMS": (20.0, "ppm"), "TOF": (40.0, "ppm"), "ITMS": (0.35, "da")}


def cfg_get(cfg, path, default=None):
    """Value at a dotted config path ("inputs.spectra_type"), or `default` if absent, null or empty."""
    node = cfg
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return default if node is None or node == "" else node


def _norm_tag(tag):
    """An Oktoberfest ``tag`` without its MSA-variant suffix ("tmtpro_msa" -> "tmtpro")."""
    t = str(tag or "").strip().lower()
    return t[:-4] if t.endswith("_msa") else t


def run_params(cfg):
    """The experiment-specific settings the report needs, resolved from the run's config.

    Resolved once and passed around, so that no panel re-derives — or hard-codes — any of them.

    :param cfg: the run's parsed config.json ({} if the run kept none)
    :return: dict with the isobaric label and its reporter-ion window (``None`` when label-free), the
        fragment match tolerance (value, unit, and where it came from), and the prediction settings
    """
    tag = _norm_tag(cfg_get(cfg, "tag", ""))
    analyzer = str(cfg_get(cfg, "inputs.instrument_type", "") or "").upper()
    tol, unit = cfg_get(cfg, "massTolerance"), str(cfg_get(cfg, "unitMassTolerance", "") or "").lower()
    if tol is not None and unit in ("ppm", "da"):
        tol, source = float(tol), "from massTolerance"
    else:  # exactly the fallback the run's own fragment annotation used
        tol, unit = ANALYZER_TOL.get(analyzer, ANALYZER_TOL["FTMS"])
        source = f"{analyzer} default" if analyzer in ANALYZER_TOL else "FTMS default, assumed"
    ion_types = cfg_get(cfg, "ion_types", CFG_DEFAULTS["ion_types"])
    return {
        "tag": tag,
        "label": LABEL_NAMES.get(tag, "none (label-free)"),
        "reporter": REPORTER_WINDOWS.get(tag),
        "tol": tol,
        "tol_unit": unit,
        "tol_text": f"{tol:g} {'ppm' if unit == 'ppm' else 'Da'}",
        "tol_source": source,
        "ion_types": "".join(str(x) for x in ion_types).lower(),
        "fragmentation": str(cfg_get(cfg, "fragmentation_method", CFG_DEFAULTS["fragmentation_method"])).upper(),
        "server": str(cfg_get(cfg, "prediction_server", CFG_DEFAULTS["prediction_server"])),
        "ssl": bool(cfg_get(cfg, "ssl", True)),
    }


def match_indices(observed_mz, predicted_mz, params):
    """Indices of `predicted_mz` within the run's fragment tolerance of one observed m/z."""
    tol = observed_mz * params["tol"] / 1e6 if params["tol_unit"] == "ppm" else params["tol"]
    return np.where(np.abs(predicted_mz - observed_mz) <= tol)[0]


# ============================== SPECTRA ==============================
GORDER = {
    "best_score": 0,
    "cutoff_accepted": 1,
    "cutoff_rejected": 2,
    "worst_score": 3,
    "decoy_high": 4,
    "decoy_low": 5,
}
B_COL, Y_COL, OTH_COL = "#1f6fb2", "#d1352b", "#8a8a85"  # b=blue, y=red, other=grey (mirror-plot convention)
UNM_COL = "#141414"  # unmatched observed peaks: near-black (was light grey) so they read clearly
PROTON = 1.007276


def select_spectra(data, n):
    """Pick which PSMs to show, n per group.

    The groups are the extremes of the score distribution, both sides of the 1% cutoff, and decoys.

    The selection runs on numpy indices, so only the ~6n chosen rows are ever materialised as a
    DataFrame. On a bulk run, building a frame of every PSM here (with its peptide and file strings)
    would cost more memory than the whole rest of the report.

    :raises ValueError: if rescore.tab lacks a column the spectra viewer cannot work without, or if
        no PSM is showable at all
    :return: (one row per distinct spectrum, {(raw, scan): [group names]})
    """
    import pandas as pd

    missing = [c for c in ("filename", "ScanNr", "Peptide", "collision_energy_aligned", "Mass") if c not in data]
    if missing:
        raise ValueError(f"rescore.tab has no {', '.join(missing)} column(s)")
    label, score, q, charge = data["label"], data["score"], data["q"], data["charge"]
    # a spectrum is only showable if we can score it, place its precursor, and re-query Koina for it
    ok = np.isfinite(score) & (charge > 0) & np.isfinite(data["collision_energy_aligned"])
    tgt = np.where(ok & (label == 1))[0]
    dec = np.where(ok & (label == -1))[0]
    conf, rej = tgt[q[tgt] <= FDR], tgt[q[tgt] > FDR]

    def top(idx, largest):
        """The n rows of `idx` with the largest (or smallest) score."""
        if idx.size == 0:
            return idx
        k = min(n, idx.size)
        v = -score[idx] if largest else score[idx]
        part = np.argpartition(v, k - 1)[:k]
        return idx[part[np.argsort(v[part])]]

    groups = {
        "best_score": top(tgt, True),
        "worst_score": top(tgt, False),
        "cutoff_accepted": top(conf, False),
        "cutoff_rejected": top(rej, True),
        "decoy_high": top(dec, True),
        "decoy_low": top(dec, False),
    }
    membership = {}
    for name, idx in groups.items():
        for row in idx.tolist():
            membership.setdefault(row, []).append(name)
    if not membership:
        raise ValueError("no scored PSMs to show")
    rows = np.array(sorted(membership), dtype=np.int64)
    ch = charge[rows].astype(int)
    mass = data["Mass"][rows]
    df = pd.DataFrame(
        {
            "filename": data["filename_levels"][data["filename"][rows]],
            "scan": data["ScanNr"][rows].astype(np.int64),
            "pep": [strip_pep(x) for x in data["Peptide_levels"][data["Peptide"][rows]]],
            "charge": ch,
            "ce": np.rint(data["collision_energy_aligned"][rows] * 100).astype(int),
            "RT": data["RT"][rows] if "RT" in data else np.full(rows.size, np.nan),
            "score": score[rows],
            "q": q[rows],
            "label": label[rows],
            "spectral_angle": data["spectral_angle"][rows] if "spectral_angle" in data else np.full(rows.size, np.nan),
            # Mass = theoretical neutral peptide mass (ExpMass in this tab is a bogus sequential index)
            "mz": (mass + ch * PROTON) / ch,
        }
    )
    df["group"] = [membership[int(r)] for r in rows]
    df["key"] = list(zip(df.filename, df.scan, strict=False))
    grp_of = {}
    for key, names in zip(df.key, df.group, strict=False):
        grp_of.setdefault(key, set()).update(names)
    grp_of = {k: sorted(v) for k, v in grp_of.items()}
    uniq = df.drop_duplicates("key").reset_index(drop=True)
    return uniq, grp_of


def _ion_type(a):
    return "b" if a and a[0] == "b" else ("y" if a and a[0] == "y" else "o")


def _stems(x, y):
    xs, ys = [], []
    for xi, yi in zip(x, y, strict=False):
        xs += [xi, xi, None]
        ys += [0, yi, None]
    return xs, ys


def _scan_of(spectrum_id):
    """Scan number encoded in an mzML spectrum id ("... scan=1234"), or None if it carries none."""
    m = re.search(r"scan=(\d+)", str(spectrum_id))
    return int(m.group(1)) if m else None


def _peaks(spectrum):
    mz, it = spectrum.get("m/z array"), spectrum.get("intensity array")
    if mz is None or it is None or len(mz) == 0:
        return None
    return np.asarray(mz, dtype=float), np.asarray(it, dtype=float)


def _read_scans(mzml_path, scans):
    """{scan: (m/z array, intensity array)} for the wanted scans of one mzML file.

    Uses the file's own index when it has one: on a bulk-size mzML, seeking to a handful of scans is
    the difference between milliseconds and a full pass over gigabytes. Falls back to reading through
    the file (stopping as soon as everything wanted has been seen) when it has no usable index.
    """
    from pyteomics import mzml as pymzml

    found = {}
    try:
        with pymzml.MzML(str(mzml_path), use_index=True) as rd:
            wanted = {}
            for spectrum_id in rd.index:
                scan = _scan_of(spectrum_id)
                if scan in scans:
                    wanted[scan] = spectrum_id
            for scan, spectrum_id in wanted.items():
                peaks = _peaks(rd.get_by_id(spectrum_id))
                if peaks is not None:
                    found[scan] = peaks
    except Exception:  # no index, an unreadable one, or an older pyteomics: fall through
        found = {}
    if len(found) == len(scans):
        return found
    with pymzml.MzML(str(mzml_path)) as rd:
        for sp in rd:
            if sp.get("ms level") != 2:
                continue
            scan = _scan_of(sp.get("id", ""))
            if scan is None:
                scan = sp.get("index", -1) + 1
            if scan in scans and scan not in found:
                peaks = _peaks(sp)
                if peaks is not None:
                    found[scan] = peaks
                if len(found) == len(scans):
                    break
    return found


def build_spectra_fig(data, paths, model, n, params):
    """Build the interactive observed-vs-predicted spectra viewer, or None if it cannot be built."""
    try:
        import pandas as pd
        from pyteomics import mzml as pymzml  # noqa: F401  # availability check for _read_scans
        from koinapy import Koina
        import plotly.graph_objects as go
    except Exception as e:
        return None, (
            f"spectra section skipped (missing dependency: {e}) — install the optional extra with "
            "`pip install oktoberfest[report]`"
        )
    if paths["spectra"] is None:
        return None, "spectra section skipped (no spectra/ dir with mzML)"
    uniq, grp_of = select_spectra(data, n)
    # koina
    kin = pd.DataFrame(
        {
            "peptide_sequences": uniq.pep.values,
            "precursor_charges": uniq.charge.values,
            "collision_energies": uniq.ce.values,
            "fragmentation_types": [params["fragmentation"]] * len(uniq),
        }
    )
    try:
        pred = Koina(model_name=model, server_url=params["server"], ssl=params["ssl"]).predict(kin)
    except Exception as e:
        return None, f"spectra section skipped (Koina query failed: {e})"

    def _dec(a):
        return a.decode() if isinstance(a, bytes) else str(a)

    lut = {}
    for key, g in pred.groupby(["peptide_sequences", "precursor_charges", "collision_energies"]):
        lut[(str(key[0]), int(key[1]), int(round(float(key[2]))))] = (
            g["mz"].to_numpy(float),
            g["intensities"].to_numpy(float),
            [_dec(a) for a in g["annotation"]],
        )

    def pred_row(r):
        return lut.get((str(r.pep), int(r.charge), int(r.ce)), (np.array([]), np.array([]), []))

    # read mzML, most-represented files first and CAPPED: the selected PSMs of a bulk run can span
    # hundreds of raw files, and every file opened is at best a seek into a multi-gigabyte mzML
    obs_by_key, opened, skipped_files, unreadable = {}, 0, 0, []
    for raw, grp in sorted(uniq.groupby("filename"), key=lambda kv: -len(kv[1])):
        mzf = paths["spectra"] / f"{raw}.mzML"
        if not mzf.exists():
            continue
        if opened >= MAX_SPECTRA_FILES:
            skipped_files += 1
            continue
        opened += 1
        try:
            scans = _read_scans(mzf, {int(sc) for sc in grp.scan})
        except Exception as e:  # a truncated or unreadable mzML costs its own spectra, not the section
            unreadable.append(f"{raw} ({e})")
            continue
        for scan, peaks in scans.items():
            obs_by_key[(raw, scan)] = peaks

    spectra = []
    for _, r in uniq.iterrows():
        obs = obs_by_key.get((r.filename, int(r.scan)))
        if obs is None:
            continue
        pmz, pint, pann = pred_row(r)
        spectra.append(
            {"row": r, "obs": obs, "pmz": pmz, "pint": pint, "pann": pann, "groups": grp_of[(r.filename, int(r.scan))]}
        )
    if not spectra:
        return None, "spectra section skipped (no observed spectra matched)"

    def gkey(d):
        return (min(GORDER.get(g, 9) for g in d["groups"]), -float(d["row"].score))

    spectra.sort(key=gkey)

    fig = go.Figure()
    spans, annos, xr = [], [], []
    for d in spectra:
        start = len(fig.data)
        r = d["row"]
        omz = np.array(d["obs"][0])
        oin = np.array(d["obs"][1], float)
        if params["reporter"]:  # label fragments, not peptide fragments -- and only in a labelled run
            lo, hi = params["reporter"]
            keep_obs = (omz < lo) | (omz > hi)
            omz, oin = omz[keep_obs], oin[keep_obs]
        pmz = np.array(d["pmz"])
        pin = np.array(d["pint"], float)
        pann = d["pann"]
        keep = (pmz > 0) & (pin > 0)
        pmz, pin, pann = pmz[keep], pin[keep], [a for a, k in zip(pann, keep, strict=False) if k]
        oin = 100 * oin / oin.max() if oin.size and oin.max() > 0 else oin
        pin = 100 * pin / pin.max() if pin.size and pin.max() > 0 else pin
        otype = np.array(["u"] * len(omz), dtype=object)
        for i, m in enumerate(omz):
            j = match_indices(m, pmz, params)
            if len(j):
                otype[i] = _ion_type(pann[j[np.argmin(np.abs(pmz[j] - m))]])
        vis0 = len(spans) == 0
        for typ, col, wid in [("u", UNM_COL, 1.0), ("b", B_COL, 1.2), ("y", Y_COL, 1.2), ("o", OTH_COL, 1.2)]:
            sel = otype == typ
            if not sel.any():
                continue
            xs, ys = _stems(omz[sel], oin[sel])
            hov = [
                f"m/z {m:.4f}<br>{v:.1f}%{' · unmatched' if typ == 'u' else ''}"
                for m, v in zip(omz[sel], oin[sel], strict=False)
            ]
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines",
                    line={"color": col, "width": wid},
                    opacity=0.55 if typ == "u" else 1.0,
                    hoverinfo="text",
                    text=sum([[h, h, ""] for h in hov], []),
                    visible=vis0,
                    showlegend=False,
                )
            )
        for typ, col in [("b", B_COL), ("y", Y_COL), ("o", OTH_COL)]:
            idx = [i for i, a in enumerate(pann) if _ion_type(a) == typ]
            if not idx:
                continue
            xs, ys = _stems(pmz[idx], -pin[idx])
            hov = [f"{pann[i]}<br>m/z {pmz[i]:.4f}<br>{pin[i]:.1f}%" for i in idx]
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines",
                    line={"color": col, "width": 1.4},
                    hoverinfo="text",
                    text=sum([[h, h, ""] for h in hov], []),
                    visible=vis0,
                    showlegend=False,
                )
            )
        spans.append((start, len(fig.data)))
        gs = " / ".join(d["groups"])
        qtxt = "decoy" if int(r.label) == -1 else f"q={float(r.q):.4g}"
        head = (
            f"<b>{r.pep}</b>  ({gs})<br>"
            f"<span style='font-size:12px'>{r.filename} · scan {int(r.scan)} · z{int(r.charge)} · "
            f"m/z {float(r.mz):.4f} · RT {float(r.RT):.1f} · NCE {int(r.ce)} · "
            f"SA {float(r.spectral_angle):.3f} · score {float(r.score):+.3f} · {qtxt}</span>"
        )
        an = [
            {
                "x": 0.0,
                "y": 1.10,
                "xref": "paper",
                "yref": "paper",
                "showarrow": False,
                "xanchor": "left",
                "text": head,
                "font": {"size": 13, "color": INK},
            }
        ]
        for i in np.argsort(pin)[::-1]:
            if pin.size == 0 or pin[i] < 8:
                break
            an.append(
                {
                    "x": float(pmz[i]),
                    "y": float(-pin[i] - 4),
                    "showarrow": False,
                    "text": pann[i],
                    "font": {
                        "size": 9,
                        "color": B_COL
                        if _ion_type(pann[i]) == "b"
                        else (Y_COL if _ion_type(pann[i]) == "y" else OTH_COL),
                    },
                }
            )
        annos.append(an)
        allmz = np.concatenate([omz, pmz]) if pmz.size else omz
        pad = (allmz.max() - allmz.min()) * 0.03 + 5
        xr.append([float(allmz.min() - pad), float(allmz.max() + pad)])

    buttons = []
    for si, d in enumerate(spectra):
        vis = [False] * len(fig.data)
        for tt in range(*spans[si]):
            vis[tt] = True
        r = d["row"]
        gs = "/".join(d["groups"])
        lab = f"{gs} | {r.pep[:22]} (score {float(r.score):+.2f}, SA {float(r.spectral_angle):.2f})"
        buttons.append(
            {
                "label": lab,
                "method": "update",
                "args": [{"visible": vis}, {"annotations": annos[si], "xaxis.range": xr[si]}],
            }
        )
    other = "".join(c for c in params["ion_types"] if c not in "by")
    legend_traces = [(B_COL, "b ion"), (Y_COL, "y ion")]
    if other:  # ion_types is configurable: a/c/x/z ions are drawn grey and have to be named
        legend_traces.append((OTH_COL, f"other ion ({'/'.join(other)})"))
    legend_traces.append((UNM_COL, "observed, unmatched"))
    for col, nm in legend_traces:
        fig.add_trace(
            go.Scatter(x=[None], y=[None], mode="lines", line={"color": col, "width": 3}, name=nm, showlegend=True)
        )
    fig.update_layout(
        updatemenus=[
            {
                "buttons": buttons,
                "x": 0.0,
                "xanchor": "left",
                "y": 1.20,
                "yanchor": "top",
                "showactive": True,
                "direction": "down",
                "font": {"size": 11},
                "bgcolor": "#f4f4f0",
                "bordercolor": "#cccccc",
            }
        ],
        annotations=annos[0],
        template="plotly_white",
        font={"family": "Inter, Helvetica, Arial, sans-serif", "size": 13, "color": INK},
        width=1050,
        height=580,
        margin={"t": 150, "b": 55, "l": 70, "r": 40},
        xaxis={"title": "m/z", "range": xr[0], "gridcolor": "#f2f2ee"},
        yaxis={
            "title": "observed ↑ / predicted ↓  (% base peak)",
            "range": [-118, 118],
            "zeroline": False,
            "gridcolor": "#eeeeea",
        },
        legend={"orientation": "h", "y": 1.16, "x": 1.0, "xanchor": "right", "font": {"size": 11}},
        hovermode="closest",
    )
    fig.add_hline(y=0, line={"color": "#cccccc", "width": 0.8})
    html = fig.to_html(full_html=False, include_plotlyjs="inline", div_id="spectra_div")
    note = (
        f"{len(spectra)} spectra across groups: "
        + ", ".join(f"{k}={sum(1 for d in spectra if k in d['groups'])}" for k in GORDER)
        + f" · matched at {params['tol_text']} ({params['tol_source']})"
    )
    if params["reporter"]:
        note += f" · {params['label']} reporter region hidden"
    if skipped_files:
        note += (
            f" · {skipped_files} further raw file(s) left unopened (cap of {MAX_SPECTRA_FILES} mzML files per report)"
        )
    if unreadable:
        note += " · unreadable mzML: " + ", ".join(unreadable[:3])
    return html, note


# ============================== native SVGs ==============================
SUMMARY_ORDER = [
    "psm_1%_FDR",
    "peptide_1%_FDR",
    "rescore_original_joint_plot_psm",
    "rescore_original_joint_plot_peptide",
    "target_vs_decoys_sa_distribution",
    "rescore_target_vs_decoys_psm_bins",
    "original_target_vs_decoys_psm_bins",
    "rescore_target_vs_decoys_peptide_bins",
    "original_target_vs_decoys_peptide_bins",
]


def collect_svgs(paths, max_gallery_files=0, max_embedded_mb=0):
    """Native Oktoberfest SVGs to embed: (summary plots, {raw file: {plot: path}}, note).

    Everything here is inlined as base64 into ONE html file, so the selection is budgeted. A bulk run
    emits three SVGs per raw file and its pooled scatter plots hold millions of vector points; without
    a cap the report would be a multi-gigabyte page that no browser opens. The summary plots are
    served first (they carry the run-level message), the per-raw gallery gets what budget is left.
    """
    svgs = {}
    for d in {paths["fdr_dir"], paths["results_dir"]}:
        if d and d.exists():
            for f in d.glob("*.svg"):
                svgs.setdefault(f.name, f)
    summary, per_raw = [], {}
    raw_re = re.compile(r"^(?P<raw>.+?)_(irt_vs_pred_rt|violin_spectral_angle_ce|mean_spectral_angle_ce)\.svg$")
    for name, path in sorted(svgs.items()):
        stem = name[:-4]
        m = raw_re.match(name)
        if m:
            per_raw.setdefault(m.group("raw"), {})[m.group(2)] = path
        elif stem in SUMMARY_ORDER or "target_vs_decoys" in stem or "FDR" in stem or "joint" in stem:
            summary.append((stem, path))
    summary.sort(key=lambda kv: SUMMARY_ORDER.index(kv[0]) if kv[0] in SUMMARY_ORDER else 99)

    notes = []
    budget = max_embedded_mb * 1e6 * 0.75 if max_embedded_mb > 0 else float("inf")  # base64 costs 4/3
    per_file_cap = MAX_SVG_MB * 1e6

    def size(path):
        try:
            return path.stat().st_size
        except OSError:
            return per_file_cap + 1  # unreadable: treat as oversized and skip it

    kept_summary, too_big, no_budget = [], 0, 0
    for stem, path in summary:
        n_bytes = size(path)
        if n_bytes > per_file_cap:
            too_big += 1
            continue
        if n_bytes > budget:
            no_budget += 1
            continue
        budget -= n_bytes
        kept_summary.append((stem, path))
    if too_big:
        notes.append(f"{too_big} native plot(s) omitted, each larger than {MAX_SVG_MB:.0f} MB")
    if no_budget:
        notes.append(f"{no_budget} further native plot(s) omitted: the {max_embedded_mb} MB embed budget was spent")

    kept_raw, n_raw = {}, len(per_raw)
    for raw in sorted(per_raw):
        if max_gallery_files > 0 and len(kept_raw) >= max_gallery_files:
            break
        plots = {k: v for k, v in per_raw[raw].items() if size(v) <= per_file_cap}
        need = sum(size(v) for v in plots.values())
        if not plots or need > budget:
            continue
        budget -= need
        kept_raw[raw] = plots
    if len(kept_raw) < n_raw:
        notes.append(
            f"gallery shows {len(kept_raw)} of {n_raw} raw files "
            f"(cap {max_gallery_files or 'none'}, {max_embedded_mb or 'unlimited'} MB embed budget)"
        )
    return kept_summary, kept_raw, "; ".join(notes)


# ============================== HTML assembly ==============================
CSS = """
:root{--ink:#161616;--ink2:#4a4a46;--muted:#8a8984;--line:#e2e0d8;--bg:#f4f2ea;--card:#fff;
  --accent:#e07b00;--accent2:#141414;--gold:#e0a82e;--green:#0a7d43;--orange:#d9691f;--red:#cf3a2e;--chip:#fbeede}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;font-family:Inter,-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
  color:var(--ink);background:var(--bg);font-size:15px;line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1240px;margin:0 auto;padding:0 24px 90px}
.masthead{background:#343131;color:#f7f2df;border:1px solid #2a2727;border-top:5px solid var(--accent);padding:26px 32px;margin:26px 0 0}
.logo{height:52px;display:block;margin:0 0 16px}
.kicker{font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:var(--gold);font-weight:700}
.masthead h1{font-size:30px;line-height:1.16;margin:6px 0 2px;font-weight:800;letter-spacing:-.015em;color:#fff}
.masthead .run{font-size:15.5px;color:#cfcabb;margin:0 0 15px}
.badges{display:flex;flex-wrap:wrap;gap:8px}
.badge{font-size:12.5px;background:#423d3d;color:#f0c98a;border:1px solid #4e4949;padding:4px 11px;font-weight:600}
.badge.grey{background:#3b3737;color:#cfcabb;border-color:#4e4949}
nav.toc{position:sticky;top:0;z-index:20;background:rgba(255,255,255,.96);backdrop-filter:blur(6px);
  border:1px solid var(--line);border-top:2px solid var(--accent2);margin:14px 0 10px;padding:6px 6px;display:flex;flex-wrap:wrap;gap:1px}
nav.toc a{font-size:13px;color:var(--ink2);text-decoration:none;padding:5px 12px}
nav.toc a:hover{background:var(--chip);color:var(--accent)}
details.sec{margin:22px 0 0;scroll-margin-top:52px}
details.sec>summary{list-style:none;cursor:pointer;display:flex;align-items:center;gap:12px;padding:11px 14px;
  background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent)}
details.sec>summary:hover{background:#fffaf2}
details.sec>summary::-webkit-details-marker{display:none}
details.sec>summary::after{content:'▸ expand';margin-left:auto;color:var(--accent);font-size:12px;font-weight:700;white-space:nowrap}
details.sec[open]>summary::after{content:'▾ collapse'}
.sectitle{font-size:20px;font-weight:800;letter-spacing:-.01em}
.secdesc-inline{color:var(--ink2);font-size:13px;font-weight:400}
.secnum{font-size:13px;font-weight:800;color:#fff;background:var(--accent2);
  min-width:27px;height:27px;display:inline-flex;align-items:center;justify-content:center;flex:none}
.secbody{padding:14px 2px 6px}
.sechead{display:flex;align-items:center;gap:12px;margin:38px 0 2px}
.sechead h2{font-size:21px;margin:0;font-weight:800;letter-spacing:-.01em}
.secdesc{color:var(--ink2);font-size:13.5px;margin:2px 0 16px 39px}
details.raw{margin:8px 0}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:10px;margin:6px 0 14px}
.kpi{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--accent);padding:14px 16px}
.kpi .v{font-size:25px;font-weight:800;letter-spacing:-.02em;line-height:1.1;color:var(--ink)}
.kpi .l{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);margin-top:7px;font-weight:700}
.kpi .s{font-size:12px;color:var(--ink2);margin-top:3px}
.kpi.pos{border-left-color:#e58a00}.kpi.pos .v{color:#c26a00}
.kpi.neg{border-left-color:var(--red)}.kpi.neg .v{color:var(--red)}
.cfg{background:var(--card);border:1px solid var(--line)}
.cfg table{width:100%;border-collapse:collapse;font-size:13.5px}
.cfg td{padding:9px 18px;border-bottom:1px solid #f0efe7}.cfg tr:last-child td{border-bottom:none}
.cfg td:first-child{color:var(--ink2);width:32%;white-space:nowrap;font-weight:600}
.cfg td:last-child{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;word-break:break-word}
.cfg tr.grp td{background:#faf8f1;color:var(--muted);font-size:10.5px;text-transform:uppercase;letter-spacing:.07em;font-weight:800;padding:7px 18px;width:auto;white-space:normal;font-family:inherit}
.cfg .cl{display:block}
.cfg .ck{display:block;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;color:var(--muted);font-weight:400}
.cfg .dflt{font-family:Inter,Helvetica,Arial,sans-serif;font-size:10px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);border:1px solid var(--line);padding:1px 5px;margin-left:8px}
.cfgnote{padding:10px 18px;border-top:1px solid var(--line);font-size:12px;color:var(--muted);line-height:1.5}
.figcard{background:var(--card);border:1px solid var(--line);padding:12px;margin:2px 0}
.figcard img{width:100%;height:auto;display:block}
.figcap{margin:13px 6px 4px;font-size:13px;color:var(--ink2);line-height:1.55;border-left:3px solid var(--accent);padding:2px 0 2px 13px}
.figcap b.h{display:block;color:var(--ink);font-size:11px;text-transform:uppercase;letter-spacing:.05em;margin-bottom:4px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:12px;align-items:start}
.svgcard{background:var(--card);border:1px solid var(--line);padding:10px}
.svgcard.half{max-width:320px}
.svgcard .t{font-size:11.5px;color:var(--ink2);font-family:ui-monospace,monospace;margin-bottom:6px;word-break:break-all}
.svgcard img{width:100%;display:block}
.spectra-wrap{background:var(--card);border:1px solid var(--line);padding:14px 16px}
.legendbar{display:flex;gap:18px;flex-wrap:wrap;font-size:12.5px;color:var(--ink2);margin:2px 0 10px;align-items:center}
.legendbar b{color:var(--ink)}
.sw{display:inline-block;width:15px;height:3px;vertical-align:middle;margin-right:6px}
details{background:var(--card);border:1px solid var(--line);padding:2px 14px;margin:8px 0}
summary{cursor:pointer;font-weight:600;font-size:14px;padding:9px 2px;font-family:ui-monospace,monospace}
summary::marker{color:var(--accent)}
.small{font-size:12px;color:var(--muted)}
footer{margin-top:46px;padding-top:18px;border-top:2px solid var(--accent2);color:var(--muted);font-size:12.5px;line-height:1.6}
.warn{color:var(--red);font-size:13px}
@media print{nav.toc{display:none}.figcard,.kpi,.cfg,.masthead,.svgcard{break-inside:avoid}}
"""


def _logo_uri():
    import base64 as _b64

    here = Path(__file__).resolve()
    candidates = [here.parent / "Oktoberfest.svg"]
    # when shipped inside the oktoberfest package, the logo lives in the repo docs
    for up in (
        here.parents[2] if len(here.parents) > 2 else here.parent,
        here.parents[3] if len(here.parents) > 3 else here.parent,
    ):
        candidates.append(up / "docs" / "_static" / "img" / "Oktoberfest.svg")
    for p in candidates:
        # best effort: an unreadable candidate just means we fall through to the next one
        with contextlib.suppress(OSError):
            if p.exists():
                return "data:image/svg+xml;base64," + _b64.b64encode(p.read_bytes()).decode()
    return ""


SECTION_DESC = {
    "summary": "Headline identification counts, the Prosit contribution, and the run configuration.",
    "yield": "Identifications accepted at every FDR threshold, with vs without Prosit — the headline comparison.",
    "sa": "How the spectral angle (observed↔predicted agreement) depends on peptide and spectrum properties.",
    "headroom": "Whether confidently-rejected target PSMs still hide identifications a better model could recover.",
    "weights": "The Percolator SVM feature weights that determine each PSM's score, with vs without Prosit.",
    "perfile": "Identifications per raw file — whether the gain is uniform or driven by a few files.",
    "movement": "How Prosit shifts each PSM's score relative to the search-features-only model.",
    "calibration": "Retention-time and mass-accuracy agreement, accepted targets vs decoys.",
    "spectra": "Observed spectra vs clean Prosit predictions for representative PSMs and decoys.",
    "native": "Figures produced by Oktoberfest itself, embedded verbatim from the run output.",
    "gallery": "Per-raw-file retention-time and collision-energy diagnostics (collapsed; click to expand).",
}


def _sechead(n, sid, title):
    return (
        f"<div class='sechead' id='{sid}'><span class='secnum'>{n}</span><h2>{title}</h2></div>"
        f"<div class='secdesc'>{SECTION_DESC.get(sid, '')}</div>"
    )


def _section(n, sid, title, body, is_open=False):
    """A collapsible section: header is the <summary>, content collapses. Summary is open by default."""
    op = " open" if is_open else ""
    desc = SECTION_DESC.get(sid, "")
    return (
        f"<details class='sec'{op} id='{sid}'><summary>"
        f"<span class='secnum'>{n}</span><span class='sectitle'>{title}</span>"
        f"<span class='secdesc-inline'>{desc}</span></summary>"
        f"<div class='secbody'>{body}</div></details>"
    )


def _kpi(v, label, s="", cls=""):
    sub = f"<div class='s'>{s}</div>" if s else ""
    return f"<div class='kpi {cls}'><div class='v'>{v}</div><div class='l'>{label}</div>{sub}</div>"


def overview_stats(paths, data):
    """Headline numbers for the summary cards, counted from the FDR-method output itself."""
    sa = data.get("spectral_angle")
    m = (data["label"] == 1) & (data["q"] <= FDR)
    s = {
        "n_rows": data["n"],
        "n_t": int((data["label"] == 1).sum()),
        "n_d": int((data["label"] == -1).sum()),
        "median_sa": float(np.median(sa[m])) if (sa is not None and m.any()) else float("nan"),
        "cutoff": float(data["score"][m].min()) if m.any() else float("nan"),
    }
    fdir, method = paths["fdr_dir"], paths["method"]
    # ids are compared as sorted hash arrays: set arithmetic stays exact and costs 8 bytes per id
    empty = np.empty(0, dtype=np.uint64)
    rid = ids_at_fdr(fdir / f"rescore.{method}.psms.txt")
    s["acc_psm"] = int(rid.size)
    pep_f = fdir / f"rescore.{method}.peptides.txt"
    rpep = ids_at_fdr(pep_f, key="peptide") if pep_f.exists() else empty
    s["acc_pep"] = int(rpep.size)
    op, opep = fdir / f"original.{method}.psms.txt", fdir / f"original.{method}.peptides.txt"
    if op.exists():
        oid = ids_at_fdr(op)
        s.update(
            gain_psm=int(rid.size - oid.size),
            gained_psm=int(np.setdiff1d(rid, oid, assume_unique=True).size),
            lost_psm=int(np.setdiff1d(oid, rid, assume_unique=True).size),
        )
    if opep.exists() and rpep.size:
        oidp = ids_at_fdr(opep, key="peptide")
        s.update(
            gain_pep=int(rpep.size - oidp.size),
            gained_pep=int(np.setdiff1d(rpep, oidp, assume_unique=True).size),
            lost_pep=int(np.setdiff1d(oidp, rpep, assume_unique=True).size),
        )
    return s


def kpi_cards(s):
    """Render the headline KPI cards from the summary statistics."""
    c = [
        _kpi(f"{s['acc_psm']:,}", "Target PSMs @1% FDR", f"of {s['n_t']:,} target PSMs"),
        _kpi(f"{s['acc_pep']:,}", "Peptides @1% FDR", ""),
    ]
    if "gain_psm" in s:
        c.append(
            _kpi(
                f"{s['gain_psm']:+,}",
                "Prosit net PSMs",
                f"+{s['gained_psm']:,} / −{s['lost_psm']:,}",
                "pos" if s["gain_psm"] >= 0 else "neg",
            )
        )
    if "gain_pep" in s:
        c.append(
            _kpi(
                f"{s['gain_pep']:+,}",
                "Prosit net peptides",
                f"+{s['gained_pep']:,} / −{s['lost_pep']:,}",
                "pos" if s["gain_pep"] >= 0 else "neg",
            )
        )
    if np.isfinite(s["median_sa"]):
        c.append(_kpi(f"{s['median_sa']:.3f}", "Median SA (accepted)"))
    c.append(_kpi(f"{s['n_d']:,}", "Decoy PSMs"))
    if np.isfinite(s["cutoff"]):
        c.append(_kpi(f"{s['cutoff']:.2f}", "Score at 1% cutoff"))
    return "<div class='kpis'>" + "".join(c) + "</div>"


#: The Summary section's parameter table: groups of (dotted config key, label). The table mirrors the
#: config file the user actually wrote — same keys, same nesting — so every row can be traced back to
#: (or copied into) a config.json, rather than being a curated selection of what happened to matter
#: for one experiment. Two rows are synthesised because no single key holds them:
#: ``__tolerance`` (massTolerance + unitMassTolerance, or the mass-analyzer fallback) and ``__label``
#: (``tag``, which is stated as "none (label-free)" instead of vanishing on an unlabelled run).
#: A group whose parent key is absent is dropped entirely — a rescoring run digests no fasta. An entry
#: may carry a third element, the unit its value is expressed in.
CONFIG_GROUPS = [
    (
        "Run",
        None,
        [
            ("type", "Workflow"),
            ("inputs.search_results_type", "Search engine"),
            ("inputs.spectra_type", "Spectra format"),
            ("inputs.instrument_type", "Mass analyzer"),
        ],
    ),
    (
        "Prediction",
        None,
        [
            ("models.intensity", "Intensity model"),
            ("models.irt", "iRT model"),
            ("models.proteotypicity", "Proteotypicity model"),
            ("prediction_server", "Prediction server"),
            ("fragmentation_method", "Fragmentation method"),
            ("ce_alignment_options.ce_range", "CE alignment range"),
            ("ce_alignment_options.use_ransac_model", "CE RANSAC model"),
        ],
    ),
    (
        "Fragment annotation",
        None,
        [
            ("__tolerance", "Fragment mass tolerance"),
            ("ion_types", "Ion types"),
            ("p_window", "Precursor exclusion window", "Da"),
            ("__label", "Isobaric label"),
            ("matching_method", "Peak matching"),
            ("matching_method_params", "Peak matching params"),
        ],
    ),
    (
        "Rescoring",
        None,
        [
            ("fdr_estimation_method", "FDR estimation"),
            ("allFeatures", "All features"),
            ("add_feature_cols", "Extra feature columns"),
            ("regressionMethod", "RT regression"),
            ("quantification", "Quantification"),
        ],
    ),
    (
        "Fasta digestion",
        "fastaDigestOptions",
        [
            ("fastaDigestOptions.enzyme", "Enzyme"),
            ("fastaDigestOptions.digestion", "Digestion"),
            ("fastaDigestOptions.missedCleavages", "Missed cleavages"),
            ("fastaDigestOptions.minLength", "Min peptide length"),
            ("fastaDigestOptions.maxLength", "Max peptide length"),
            ("fastaDigestOptions.db", "Database"),
        ],
    ),
    (
        "Custom modifications",
        "inputs.custom_modifications",
        [
            ("inputs.custom_modifications.static_mods", "Static mods"),
            ("inputs.custom_modifications.var_mods", "Variable mods"),
        ],
    ),
]


def _fmt_cfg(v):
    """A config value as it reads in a config.json, flattened to one line."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, dict):
        return ", ".join(f"{k}: {_fmt_cfg(x)}" for k, x in v.items()) if v else "—"
    if isinstance(v, (list, tuple)):
        return ", ".join(_fmt_cfg(x) for x in v) if len(v) else "—"
    return str(v)


def config_rows(cfg, params):
    """The parameter table's rows: (group, label, config key, value, fallback note).

    The note is empty for a value the config actually set, and otherwise names the fallback the run
    used in its place ("default", or the mass analyzer whose tolerance stood in for a missing one).

    :param cfg: the run's parsed config.json
    :param params: resolved run settings from :py:func:`run_params`
    :return: list of row tuples, groups in order, empty groups already dropped
    """
    rows = []
    for group, parent, keys in CONFIG_GROUPS:
        if parent is not None and cfg_get(cfg, parent) is None:
            continue  # section the run never used (no fasta digestion, no custom mods)
        for key, label, *unit in keys:
            if key == "__tolerance":
                explicit = params["tol_source"].startswith("from")
                rows.append(
                    (
                        group,
                        label,
                        "massTolerance / unitMassTolerance",
                        params["tol_text"],
                        "" if explicit else params["tol_source"],
                    )
                )
                continue
            if key == "__label":
                rows.append((group, label, "tag", params["label"], "" if params["tag"] else "default"))
                continue
            given = cfg_get(cfg, key)
            if given is None and key not in CFG_DEFAULTS:
                continue  # optional and unset (e.g. no proteotypicity model, no instrument_type)
            value = given if given is not None else CFG_DEFAULTS[key]
            shown = _fmt_cfg(value) + (f" {unit[0]}" if unit and value != "" else "")
            rows.append((group, label, key, shown, "" if given is not None else "default"))
    return rows


def config_table(cfg, params):
    """The Summary section's parameter table, or a warning if the run kept no config."""
    if not cfg:
        return (
            "<div class='cfg'><div class='cfgnote warn'>No config.json was found next to this run, "
            "so its parameters cannot be shown. Everything below is read from the run output "
            "itself; the spectra viewer falls back to Oktoberfest's own defaults.</div></div>"
        )
    html, seen = [], None
    for group, label, key, value, note in config_rows(cfg, params):
        if group != seen:
            html.append(f"<tr class='grp'><td colspan='2'>{group}</td></tr>")
            seen = group
        badge = f"<span class='dflt'>{note}</span>" if note else ""
        html.append(
            f"<tr><td><span class='cl'>{label}</span><span class='ck'>{key}</span></td><td>{value}{badge}</td></tr>"
        )
    return (
        "<div class='cfg'><table>" + "".join(html) + "</table>"
        "<div class='cfgnote'>Keys are the paths in this run's <b>config.json</b>, so every row can "
        "be read straight off — or copied straight into — a config. A value carrying a "
        "<span class='dflt'>default</span> marker was <i>not</i> set there: it is what Oktoberfest "
        "fell back to, not a choice the config records.</div></div>"
    )


def build_report_html(
    paths,
    data,
    model,
    run_name,
    cfg,
    params,
    stats,
    sections,
    spectra_html,
    spectra_note,
    summary_svgs,
    per_raw,
    want_spectra,
    svg_note="",
):
    """Assemble every section into the single self-contained HTML document."""
    nav = [("summary", "Summary")] + [(sid, title) for sid, title, _, _ in sections]
    if spectra_html:
        nav.append(("spectra", "Spectra"))
    if summary_svgs:
        nav.append(("native", "Native plots"))
    if per_raw:
        nav.append(("gallery", "Per-raw gallery"))
    engine = (cfg.get("inputs", {}) or {}).get("search_results_type", "")
    labelled = [params["label"]] if params["tag"] else []
    badges = "".join(f"<span class='badge'>{b}</span>" for b in [model] + ([engine] if engine else []) + labelled)
    badges += f"<span class='badge grey'>generated {datetime.now():%Y-%m-%d %H:%M}</span>"

    h = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>Oktoberfest report — {run_name}</title><style>{CSS}</style></head><body><div class='wrap'>",
    ]
    logo = _logo_uri()
    logo_html = f"<img class='logo' src='{logo}' alt='Oktoberfest'>" if logo else ""
    h.append(
        f"<div class='masthead'>{logo_html}<div class='kicker'>Rescoring report</div>"
        f"<h1>{run_name}</h1><div class='run'>Spectral-prediction PSM rescoring — diagnostic report</div>"
        f"<div class='badges'>{badges}</div></div>"
    )
    h.append("<nav class='toc'>" + "".join(f"<a href='#{sid}'>{t}</a>" for sid, t in nav) + "</nav>")

    n = 1
    h.append(
        _section(n, "summary", "Summary", kpi_cards(stats) + config_table(cfg, params), is_open=True)
    )  # only Summary open
    for sid, title, cap, uri in sections:
        n += 1
        body = (
            f"<div class='figcard'><img src='{uri}' alt='{title}'>"
            f"<div class='figcap'><b class='h'>What this shows</b>{cap}</div></div>"
        )
        h.append(_section(n, sid, title, body))
    if spectra_html:
        n += 1
        other = "".join(c for c in params["ion_types"] if c not in "by")
        other_sw = (
            (f"<span><span class='sw' style='background:{OTH_COL}'></span>other ion ({'/'.join(other)})</span>")
            if other
            else ""
        )
        hidden = (
            (f" · {params['label']} reporter ions {params['reporter'][0]:.1f}–{params['reporter'][1]:.1f} m/z hidden")
            if params["reporter"]
            else ""
        )
        legend = (
            "<div class='legendbar'>"
            f"<span><span class='sw' style='background:{B_COL}'></span><b>b</b> ion (matched)</span>"
            f"<span><span class='sw' style='background:{Y_COL}'></span><b>y</b> ion (matched)</span>"
            f"{other_sw}"
            f"<span><span class='sw' style='background:{UNM_COL}'></span>observed, unmatched</span>"
            f"<span>observed ↑ / predicted ↓{hidden} · {params['tol_text']} match window</span></div>"
        )
        body = (
            "<div class='spectra-wrap'>"
            + legend
            + "<div class='figcap' style='margin:0 0 8px'><b class='h'>How to read</b>"
            "Observed spectrum (up, from mzML) mirrored against the clean Koina re-query (down). "
            "Use the dropdown to switch PSM. Each header carries the full info line: "
            f"raw · scan · peptide(+mods) · charge · m/z · RT · NCE · SA · score · q. Groups: best/worst {paths['method']} "
            "score, both sides of the score≈0 cutoff, plus high- and low-scoring decoys.</div>"
            f"<div class='small' style='margin:0 0 8px'>{spectra_note}</div>" + spectra_html + "</div>"
        )
        h.append(_section(n, "spectra", "Spectra viewer", body))
    elif want_spectra:
        n += 1
        h.append(
            _section(
                n, "spectra", "Spectra viewer", f"<div class='figcard'><div class='warn'>{spectra_note}</div></div>"
            )
        )
    if summary_svgs:
        n += 1
        half = {"psm_1%_FDR", "peptide_1%_FDR"}  # count plots -> render at half size
        cells = "".join(
            f"<div class='svgcard{' half' if stem in half else ''}'><div class='t'>{stem}</div>"
            f"<img src='{svg_to_uri(p)}'></div>"
            for stem, p in summary_svgs
        )
        h.append(_section(n, "native", "Native Oktoberfest plots", f"<div class='grid2'>{cells}</div>"))
    if per_raw:
        n += 1
        n_full = sum(1 for r in per_raw.values() if len(r) >= 3)
        notes = []
        if n_full < len(per_raw):
            notes.append(
                f"Oktoberfest emitted the iRT-vs-RT plot for only {n_full}/{len(per_raw)} of the raw files "
                "shown, so some entries have 2 plots (CE violin + mean) rather than 3 — the missing SVG "
                "was never produced by the run."
            )
        if svg_note:
            notes.append(svg_note[0].upper() + svg_note[1:] + ".")
        note = ("<div class='secdesc' style='margin:0 0 10px 0'>Note: " + " ".join(notes) + "</div>") if notes else ""
        body = note + "".join(
            f"<details class='raw'><summary>{raw} · {len(per_raw[raw])} plots</summary>"
            f"<div class='grid2' style='margin:8px 0'>"
            + "".join(
                f"<div class='svgcard'><div class='t'>{k}</div><img src='{svg_to_uri(v)}'></div>"
                for k, v in sorted(per_raw[raw].items())
            )
            + "</div></details>"
            for raw in sorted(per_raw)
        )
        h.append(_section(n, "gallery", "Per-raw gallery", body))
    h.append(
        f"<footer>Generated by <b>oktoberfest.plotting.investigate</b> · source: {paths['fdr_dir']} "
        f"({paths['method']}) · "
        f"{datetime.now():%Y-%m-%d %H:%M}.<br>Non-native panels define every colour, mark, unit and n in-figure. "
        "Spectra are observed mzML vs a clean Koina re-query (not the rescoring's stored predictions). "
        "A diagnostic aid alongside — not a replacement for — the primary Oktoberfest outputs.</footer>"
    )
    # open a collapsed section when its nav link (or any in-page anchor) targets it
    h.append(
        "<script>function _openHash(){var el=document.getElementById((location.hash||'').slice(1));"
        "if(el&&el.tagName==='DETAILS'){el.open=true;el.scrollIntoView();}}"
        "window.addEventListener('hashchange',_openHash);window.addEventListener('load',_openHash);</script>"
    )
    h.append("</div></body></html>")
    return "\n".join(h)


def build_report(
    out_dir,
    out_html=None,
    n_per_group=20,
    want_spectra=True,
    want_pdf=False,
    max_psms=0,
    max_gallery_files=0,
    max_embedded_mb=0,
    log=print,
):
    """Build the HTML report for an Oktoberfest run.

    :param out_dir: the run's out/ folder, or directly the percolator/mokapot output dir
    :param out_html: where to write the report; defaults to investigate_<run>.html next to the run
    :param n_per_group: spectra shown per selection group in the spectra viewer
    :param want_spectra: build the interactive spectra viewer (needs mzML files and Koina)
    :param want_pdf: additionally render a printable PDF of the report (needs headless Chrome)
    :param max_psms: refuse to build the report above this many PSMs (0 = no limit)
    :param max_gallery_files: raw files shown in the per-raw SVG gallery (0 = no limit)
    :param max_embedded_mb: budget for all embedded native SVGs together (0 = no limit)
    :param log: where progress goes; the pipeline passes its logger
    :return: the path of the written HTML report

    Individual figures and the spectra section are already guarded here; callers that must never fail
    (e.g. inside the pipeline) should still use :py:func:`build_report_safe`.
    """
    paths = find_paths(out_dir)
    cfg = {}
    model = "Prosit_2025_intensity_40PTM"
    if paths["config"]:
        try:
            cfg = json.loads(Path(paths["config"]).read_text())
            model = cfg.get("models", {}).get("intensity", model)
        except Exception:
            cfg = {}
    run_name = paths["run_dir"].name
    params = run_params(cfg)  # label, tolerance, fragmentation, server: read once, never assumed
    engine = (cfg.get("inputs", {}) or {}).get("search_results_type", "")
    RUN_CTX["suffix"] = "  ·  ".join(x for x in [engine, model, run_name] if x)  # adaptive figure subtitles
    out_html = Path(out_html) if out_html else (paths["run_dir"] / f"investigate_{run_name}.html")

    log(
        f"[investigate] run={run_name} dir={paths['fdr_dir']} fdr={paths['method']} model={model} "
        f"label={params['tag'] or 'none'} tol={params['tol_text']} ({params['tol_source']})"
    )
    data = load_data(paths, max_psms=max_psms, log=log)
    log(f"[investigate]   loaded {data['n']:,} PSMs")

    sections = []
    for fn in [fig_yield, fig_movement, fig_per_file, fig_weights, fig_sa_features, fig_headroom, fig_calibration]:
        try:
            res = fn(data, paths)
        except Exception as e:
            log(f"[investigate]   skip {fn.__name__}: {e}")
            res = None
        if res:
            sections.append(res)

    spectra_html, spectra_note = None, "spectra section disabled"
    if want_spectra:
        try:
            spectra_html, spectra_note = build_spectra_fig(data, paths, model, n_per_group, params)
        except Exception as e:
            spectra_html, spectra_note = None, f"spectra section skipped ({e})"
        log(f"[investigate]   {spectra_note}")

    summary_svgs, per_raw, svg_note = collect_svgs(paths, max_gallery_files, max_embedded_mb)
    if svg_note:
        log(f"[investigate]   {svg_note}")
    try:
        stats = overview_stats(paths, data)
    except Exception as e:  # a summary without its headline numbers still beats no report at all
        log(f"[investigate]   summary numbers unavailable: {e}")
        acc = (data["label"] == 1) & (data["q"] <= FDR)
        stats = {
            "n_rows": data["n"],
            "n_t": int((data["label"] == 1).sum()),
            "n_d": int((data["label"] == -1).sum()),
            "median_sa": float("nan"),
            "cutoff": float("nan"),
            "acc_psm": int(acc.sum()),
            "acc_pep": 0,
        }
    html = build_report_html(
        paths,
        data,
        model,
        run_name,
        cfg,
        params,
        stats,
        sections,
        spectra_html,
        spectra_note,
        summary_svgs,
        per_raw,
        want_spectra,
        svg_note,
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html)
    log(
        f"[investigate] wrote {out_html} ({out_html.stat().st_size / 1e6:.1f} MB, {len(sections)} figs, "
        f"{'spectra' if spectra_html else 'no spectra'}, {len(per_raw)} raw galleries)"
    )

    if want_pdf:
        from oktoberfest.plotting.report_pdf import html_to_pdf

        html_to_pdf(out_html, log=log)  # guarded internally: a missing Chrome costs the PDF, not the HTML
    return out_html


def build_report_safe(
    out_dir,
    out_html=None,
    n_per_group=20,
    want_spectra=True,
    want_pdf=False,
    max_psms=0,
    max_gallery_files=0,
    max_embedded_mb=0,
    log=print,
):
    """Never-raises wrapper for pipeline use: catches EVERYTHING and returns the path or None.

    Same arguments as :py:func:`build_report`.
    """
    try:
        return build_report(
            out_dir,
            out_html=out_html,
            n_per_group=n_per_group,
            want_spectra=want_spectra,
            want_pdf=want_pdf,
            max_psms=max_psms,
            max_gallery_files=max_gallery_files,
            max_embedded_mb=max_embedded_mb,
            log=log,
        )
    except Exception as e:  # noqa: BLE001 - deliberately swallow all errors so a report never breaks a run
        # the caller-supplied log callback may itself raise; that must not mask the real failure
        with contextlib.suppress(Exception):
            log(f"[investigate] report generation failed, skipped: {e}")
        return None


def main():
    """Build a report from the command line, for a run that has already finished."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir", help="an Oktoberfest run's out/ folder, or its percolator/mokapot dir")
    ap.add_argument("out_html", nargs="?", default=None)
    ap.add_argument("--n-per-group", type=int, default=20)
    ap.add_argument("--no-spectra", action="store_true", help="skip the interactive spectra viewer")
    ap.add_argument("--pdf", action="store_true", help="also render a printable PDF (needs headless Chrome)")
    ap.add_argument("--max-psms", type=int, default=0, help="refuse runs larger than this (0 = no limit)")
    ap.add_argument("--max-gallery-files", type=int, default=25, help="raw files in the per-raw gallery")
    ap.add_argument("--max-embedded-mb", type=int, default=60, help="embed budget for native SVGs")
    a = ap.parse_args()
    build_report(
        a.out_dir,
        a.out_html,
        a.n_per_group,
        want_spectra=not a.no_spectra,
        want_pdf=a.pdf,
        max_psms=a.max_psms,
        max_gallery_files=a.max_gallery_files,
        max_embedded_mb=a.max_embedded_mb,
    )


if __name__ == "__main__":
    main()
