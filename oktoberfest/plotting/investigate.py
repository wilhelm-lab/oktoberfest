#!/usr/bin/env python3
"""Oktoberfest RUN INVESTIGATOR — one script, one self-contained HTML.

Point it at an Oktoberfest Rescoring output folder (`<run>/out`, or directly a percolator dir) and it
builds a single HTML report bundling:
  * identification yield vs FDR threshold, +Prosit vs no-Prosit (the headline comparison),
  * rescore-vs-original per-PSM score movement, coloured by what was gained/lost/kept,
  * identifications per raw file — whether the gain is uniform across files,
  * Percolator feature weights (no-Prosit vs +Prosit),
  * SA diagnostics (spectral angle vs the properties that move it),
  * rescoring headroom (rejected targets vs decoys),
  * RT & mass-error calibration, accepted targets vs decoys,
  * native Oktoberfest plots (gain/loss bars, target-vs-decoy, SA distribution, rescore-vs-original
    joint) + a collapsible gallery of every per-raw SVG (iRT-vs-RT, CE spectral-angle violins),
  * an interactive spectra viewer: observed vs clean Koina prediction, score-based groups + decoys,
    unmatched observed peaks drawn near-black, full info per spectrum (raw/scan/peptide/charge/RT/CE/
    SA/score/q/m-over-z).

Generic across search engines: every optional feature column / SVG is used only if present.

CLARITY RULE (enforced on every non-native panel): a legend or annotation DEFINES every colour and
mark, every axis has units, every distribution shows n, and each panel has a one-line "what is this".
EARN-YOUR-PLACE RULE: a panel that restates a native Oktoberfest plot, a KPI card or a neighbouring
panel does not go in. Deliberate omissions are recorded as comments where they would otherwise be
re-added out of habit (see fig_sa_features, fig_calibration).

Usage (cluster compute node, mp26_test env — needs mzML + Koina for the spectra section):
    python oktoberfest_investigate.py <run>/out [out.html] [--n-per-group 20] [--no-spectra]
Runs locally too (on a synced percolator dir); the spectra section auto-skips if mzML/Koina absent.
For a printable PDF of an already-built report, see tools/report_to_pdf.py.
"""
import argparse
import base64
import glob
import io
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
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
COMMON, GAINED, LOST = "#115795", "#007D3E", "#E17224"   # native Oktoberfest gains/losses colours
TARGET_C, DECOY_C = "#2f9e44", "#8f8e88"                 # target = green, decoy = neutral grey
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#c9c9c4", "axes.linewidth": 0.9,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.labelcolor": INK2})
FDR = 0.01


def fig_to_uri(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def svg_to_uri(path):
    return "data:image/svg+xml;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def grid(ax, axis="both"):
    ax.grid(axis=axis, color=GRID, lw=0.7); ax.set_axisbelow(True)


# ============================== data loading ==============================
def find_paths(base):
    base = Path(base).resolve()
    # percolator dir = the dir that actually holds rescore.percolator.psms.txt
    cand = [base, base / "results" / "percolator", base / "percolator"]
    cand += [p.parent for p in base.rglob("rescore.percolator.psms.txt")]
    perc = next((c for c in cand if (c / "rescore.percolator.psms.txt").exists()), None)
    if perc is None:
        raise FileNotFoundError(f"could not find rescore.percolator.psms.txt under {base}")
    # spectra dir + config: look near an `out` root
    out_root = base if (base / "spectra").exists() else perc.parent.parent
    spectra = next((d for d in [out_root / "spectra", base / "spectra"] if d.exists()), None)
    run_dir = out_root.parent
    config = next((c for c in [run_dir / "config.json", out_root / "config.json",
                               base / "config.json"] if c.exists()), None)
    results = perc.parent  # out/results (holds the CE-violin SVGs) is perc.parent
    return dict(perc=perc, spectra=spectra, config=config, run_dir=run_dir,
                out_root=out_root, results_dir=out_root / "results" if (out_root / "results").exists() else results)


def read_tab(path, wanted):
    """Ragged-safe manual parse; returns dict col->np.array for the wanted cols that EXIST."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    cols = [c for c in wanted if c in header]
    idx = [header.index(c) for c in cols]
    mx = max(idx)
    buf = {c: [] for c in cols}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= mx:
                continue
            for c, i in zip(cols, idx):
                buf[c].append(p[i])
    out = {}
    for c in cols:
        if c in ("SpecId", "Peptide", "filename", "Proteins"):
            out[c] = np.array(buf[c], dtype=object)
        else:
            out[c] = np.array(buf[c], float)
    return out, header


def load_perc(path):
    """PSMId -> (score, q). Manual parse (proteinIds carries extra tabs)."""
    d = {}
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        pi, si, qi = h.index("PSMId"), h.index("score"), h.index("q-value")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d[p[pi]] = (float(p[si]), float(p[qi]))
    return d


def ids_at_fdr(path, keycol, fdr=FDR):
    ids = set()
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        ki, qi = h.index(keycol), h.index("q-value")
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                if float(p[qi]) <= fdr:
                    ids.add(p[ki])
            except (IndexError, ValueError):
                pass
    return ids


def parse_weights(path):
    """Percolator weights.csv -> (feature_names, mean normalized weight per feature).
    Layout after the 3 comment lines: (names, normalized, raw) repeated per CV bin."""
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
    return names, dict(zip(names, mean))


def strip_pep(p):
    return p[2:-2] if isinstance(p, str) and p.startswith("_.") and p.endswith("._") else p


# ============================== load everything ==============================
def load_data(P):
    perc = P["perc"]
    wanted = ["SpecId", "Label", "ScanNr", "filename", "ExpMass", "Mass", "RT", "iRT", "pred_RT",
              "abs_rt_diff", "collision_energy_aligned", "missedCleavages", "sequence_length", "KR",
              "delta_mass_ppm", "mean_ppm_error", "max_ppm_error", "log10_evalue", "pearson_corr",
              "spectral_angle", "spectral_angle_no_b1", "spectral_angle_noise_aware",
              "fraction_observed_and_predicted", "count_observed_and_predicted", "count_predicted",
              "Peptide", "Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "Charge6"]
    t, header = read_tab(perc / "rescore.tab", wanted)
    n = len(t["SpecId"])
    charge = np.zeros(n, int)
    for c in range(1, 7):
        key = f"Charge{c}"
        if key in t:
            charge[t[key] == 1] = c
    label = t["Label"].astype(int)

    tgt = load_perc(perc / "rescore.percolator.psms.txt")
    decp = perc / "rescore.percolator.decoy.psms.txt"
    dec = load_perc(decp) if decp.exists() else {}
    score = np.full(n, np.nan); q = np.full(n, np.nan)
    for i, sid in enumerate(t["SpecId"]):
        rec = tgt.get(sid) if label[i] == 1 else dec.get(sid)
        if rec is not None:
            score[i] = rec[0]; q[i] = rec[1] if label[i] == 1 else np.nan

    origp = perc / "original.percolator.psms.txt"
    origd = load_perc(origp) if origp.exists() else {}
    origdec = perc / "original.percolator.decoy.psms.txt"
    origdd = load_perc(origdec) if origdec.exists() else {}
    oscore = np.full(n, np.nan); oq = np.full(n, np.nan)
    for i, sid in enumerate(t["SpecId"]):
        rec = origd.get(sid) if label[i] == 1 else origdd.get(sid)
        if rec is not None:
            oscore[i] = rec[0]
            if label[i] == 1:
                oq[i] = rec[1]

    t.update(dict(charge=charge, label=label, score=score, q=q, oscore=oscore, oq=oq, header=header, n=n))
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
    """Native-Oktoberfest-style scatter: every point plotted, coloured by local 2-D density
    (warm cmap, log scale, dense drawn on top) so sparse vs dense regions are legible."""
    from scipy.stats import spearmanr
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    if x.size < 20:
        ax.axis("off"); return
    lo, hi = np.percentile(x, clip); k = (x >= lo) & (x <= hi); x, y = x[k], y[k]
    n_true = x.size
    nb = 64
    H, xe, ye = np.histogram2d(x, y, bins=nb)
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, nb - 1)
    iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, nb - 1)
    dens = H[ix, iy]
    # trend + correlation on the FULL data
    r, _ = spearmanr(x, y)
    mx, my = [], []
    if trend:
        bins = np.linspace(x.min(), x.max(), 12); bi = np.clip(np.digitize(x, bins), 1, len(bins) - 1)
        for b in range(1, len(bins)):
            s = bi == b
            if s.sum() > 20:
                mx.append(x[s].mean()); my.append(np.median(y[s]))
    # scatter (subsample for file size, keeping the density ordering so dense stays on top)
    o = np.argsort(dens, kind="stable"); xs, ys, ds = x[o], y[o], dens[o]
    if xs.size > 55000:
        sel = np.linspace(0, xs.size - 1, 55000).astype(int); xs, ys, ds = xs[sel], ys[sel], ds[sel]
    sc = ax.scatter(xs, ys, c=np.log1p(ds), cmap="YlOrRd", s=5, linewidths=0, rasterized=True)
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_ticks([]); cb.ax.set_ylabel("point density", fontsize=8, color=INK2)
    cb.ax.text(0.5, 1.01, "dense", transform=cb.ax.transAxes, ha="center", va="bottom", fontsize=7, color=INK2)
    cb.ax.text(0.5, -0.01, "sparse", transform=cb.ax.transAxes, ha="center", va="top", fontsize=7, color=INK2)
    if trend:
        ax.plot(mx, my, "-o", color=INK, lw=1.8, ms=3, label="binned median")
        ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.text(0.97, 0.06, f"Spearman r={r:+.2f}\nn={n_true:,}", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color=INK, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))
    _caption(ax, title)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylab); ax.set_xlabel(xlabel)
    grid(ax)


def _violins(ax, groups, positions, color, width=0.75):
    v = ax.violinplot(groups, positions=positions, widths=width, showextrema=False)
    for b in v["bodies"]:
        b.set_facecolor(color); b.set_edgecolor(INK2); b.set_alpha(0.55); b.set_linewidth(1)
    for g, p in zip(groups, positions):
        med = np.median(g)
        ax.hlines(med, p - width / 3, p + width / 3, color=INK, lw=1.8)
        ax.text(p, 1.03, f"{med:.2f}", ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")


def _split_violin(ax, left, right, colL, colR, pos=0.0, width=0.9):
    """One violin split in two: left half = `left` distribution, right half = `right`."""
    for data, side, col in [(left, "left", colL), (right, "right", colR)]:
        v = ax.violinplot([data], positions=[pos], widths=width, showextrema=False)
        b = v["bodies"][0]
        verts = b.get_paths()[0].vertices
        if side == "left":
            verts[:, 0] = np.clip(verts[:, 0], -np.inf, pos)
        else:
            verts[:, 0] = np.clip(verts[:, 0], pos, np.inf)
        b.set_facecolor(col); b.set_edgecolor(INK2); b.set_alpha(0.6); b.set_linewidth(1)
    ax.hlines(np.median(left), pos - width / 2, pos, color=INK, lw=2)
    ax.hlines(np.median(right), pos, pos + width / 2, color=INK, lw=2)
    ax.text(pos - width / 4, 1.03, f"{np.median(left):.2f}", ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")
    ax.text(pos + width / 4, 1.03, f"{np.median(right):.2f}", ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")


def fig_sa_features(D, P):
    acc = (D["label"] == 1) & (D["q"] <= FDR)
    sa = D["spectral_angle"]; saa = sa[acc]
    ch = D["charge"][acc]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5)); A = axes.ravel()
    # A0: SA by precursor charge — split violin (2+ | 3+) when both present, else side-by-side violins
    c2, c3 = saa[ch == 2], saa[ch == 3]
    if len(c2) > 30 and len(c3) > 30:
        _split_violin(A[0], c2, c3, BLUE, ORANGE)
        A[0].set_xticks([]); A[0].set_xlim(-0.7, 0.7)
        A[0].legend(handles=[Patch(facecolor=BLUE, alpha=0.6, edgecolor=INK2, label=f"charge 2+  (n={len(c2):,})"),
                             Patch(facecolor=ORANGE, alpha=0.6, edgecolor=INK2, label=f"charge 3+  (n={len(c3):,})")],
                    frameon=False, fontsize=9, loc="lower center")
        _caption(A[0], "SA by precursor charge  (split violin: 2+ | 3+)")
    else:
        chs = [c for c in [1, 2, 3, 4, 5, 6] if (ch == c).sum() > 30]
        _violins(A[0], [saa[ch == c] for c in chs], list(range(len(chs))), BLUE)
        A[0].set_xticks(range(len(chs))); A[0].set_xticklabels([f"{c}+\n(n={int((ch == c).sum()):,})" for c in chs])
        A[0].set_xlabel("precursor charge"); _caption(A[0], "SA by precursor charge  (violins)")
    A[0].set_ylim(0, 1.12); A[0].set_ylabel("spectral angle"); grid(A[0], "y")
    # A1: SA by missed cleavages — violins
    mc = np.clip(D["missedCleavages"][acc], 0, 2).astype(int)
    mcg = [saa[mc == k] for k in [0, 1, 2]]
    _violins(A[1], mcg, [0, 1, 2], GOLD)
    A[1].set_xticks([0, 1, 2]); A[1].set_xticklabels([f"{'2+' if k == 2 else k}\n(n={len(mcg[k]):,})" for k in [0, 1, 2]])
    A[1].set_ylim(0, 1.12); A[1].set_ylabel("spectral angle"); A[1].set_xlabel("missed cleavages")
    _caption(A[1], "SA by missed cleavages  (violins)"); grid(A[1], "y")
    # A2..A3: the two continuous covariates that actually move SA. Deliberately NOT plotted, because each
    # was flat or duplicated a panel we keep: retention time (Spearman ~+0.05 — no effect), matched-fragment
    # COUNT and peptide LENGTH (collinear with fragment coverage / precursor mass but weaker), and SA-vs-RT-
    # residual (weak cross-term; the RT residual gets its own target-vs-decoy panel under Calibration).
    panels = [("fraction_observed_and_predicted", "fragment coverage", "fraction obs & pred"),
              ("Mass", "precursor mass", "precursor mass (Da)")]
    ai = 2
    for col, title, lab in panels:
        if ai >= 4:
            break
        if col not in D or not np.isfinite(D[col][acc]).any():
            continue
        density_scatter(A[ai], D[col][acc], saa, title, lab); ai += 1
    for j in range(ai, 4):
        A[j].axis("off")
    _sup(fig, "SA diagnostics — spectral angle vs features (confident targets, q≤1%)")
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=3.2, w_pad=2.2)
    cap = ("Confident target PSMs (q≤1%), showing the four properties that measurably move the spectral angle. "
           "Violin panels: shaded = SA distribution, black bar = median (value above); SA-by-charge is a split "
           "violin — left half 2+, right half 3+. Scatter panels plot every PSM coloured by local point density "
           "(warm colourbar = dense), with a black binned-median line; Spearman r and n annotated, x clipped to "
           "[p1,p99]. Reading: higher charge, missed cleavages and heavier precursors all lower SA, while fragment "
           "coverage raises it — i.e. SA degrades exactly where the fragmentation model is least constrained.")
    return "sa", "SA diagnostics", cap, fig_to_uri(fig)


def fig_headroom(D, P):
    sa, label, score, q = D["spectral_angle"], D["label"], D["score"], D["q"]
    has = np.isfinite(score)
    acc = (label == 1) & (q <= FDR); rej = (label == 1) & (q > FDR); isd = (label == -1)
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8)); bins = np.linspace(0, 1, 60)
    # (a)
    ax = axes[0]
    for m, col, lbl in [(acc, GREEN, f"accepted targets (q≤1%), n={int(acc.sum()):,}"),
                        (rej, ORANGE, f"rejected targets (q>1%), n={int(rej.sum()):,}"),
                        (isd, INK2, f"decoys, n={int(isd.sum()):,}")]:
        ax.hist(sa[m], bins=bins, density=True, histtype="step", lw=2, color=col, label=lbl)
    ax.set_xlabel("spectral angle"); ax.set_ylabel("density (area=1)")
    ax.legend(frameon=False, fontsize=9); _caption(ax, "(a) SA: accepted vs rejected targets vs decoys"); grid(ax)
    # (b) — the quantitative version of (a): would an SA cut recover anything the decoys don't also pass?
    ax = axes[1]
    thr = [0.5, 0.6, 0.7, 0.8]
    rc = [int((rej & (sa > th)).sum()) for th in thr]; dc = [int((isd & (sa > th)).sum()) for th in thr]
    x = np.arange(len(thr)); w = 0.38
    ax.bar(x - w / 2, rc, w, color=ORANGE, label="rejected targets (q>1%)")
    ax.bar(x + w / 2, dc, w, color=INK2, label="decoys")
    for i, (r, d) in enumerate(zip(rc, dc)):
        ax.text(i, max(r, d) + max(rc) * 0.02, f"net {r - d:+,}", ha="center", va="bottom",
                fontsize=10, color=GREEN if r - d > 0 else RED, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([f"SA>{th}" for th in thr]); ax.set_ylabel("PSM count")
    ax.set_ylim(0, max(rc) * 1.2); ax.legend(frameon=False, fontsize=9)
    _caption(ax, "(b) high-SA rejected targets vs high-SA decoys"); grid(ax, "y")
    _sup(fig, "Rescoring headroom — target PSMs NOT accepted at 1% FDR", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    # The caption states how to READ the net numbers, not what they come out as: the sign is run-dependent
    # (it is positive on some search-engine/model combinations and negative on others), so asserting a
    # direction here would make the report wrong on half the runs it is generated for.
    cap = ("Do we leave good identifications behind? (a) the SA distribution of REJECTED targets, compared against "
           "the DECOYS: where the two overlap, the rejects are indistinguishable from noise; any excess of rejected "
           "targets over decoys at high SA is the part that could still be real. (b) the same comparison in counts. "
           "Decoys estimate how many false PSMs an SA cut would also let through, so net = (rejected − decoy) "
           "approximates the true targets such a cut could recover: a large positive net means real headroom that "
           "a better spectral model could reach, while a net near zero or negative means an SA cut would buy nothing "
           "the FDR has not already priced in.")
    return "headroom", "Rescoring headroom", cap, fig_to_uri(fig)


def fig_weights(D, P):
    resc = P["perc"] / "rescore.percolator.weights.csv"
    orig = P["perc"] / "original.percolator.weights.csv"
    if not resc.exists():
        return None
    rn, rw = parse_weights(resc)
    on, ow = parse_weights(orig) if orig.exists() else ([], {})
    two = bool(ow)
    fig, axes = plt.subplots(1, 2 if two else 1, figsize=(14.5 if two else 8, 9.2), squeeze=False)
    def panel(ax, w, title, added=set()):
        items = [(k, v) for k, v in w.items() if abs(v) > 1e-4]  # drop exact-zero (unused) features
        items = sorted(items, key=lambda kv: abs(kv[1]), reverse=True)[:22][::-1]
        names = [k for k, _ in items]; vals = [v for _, v in items]
        cols = [(GREEN if names[i] in added else BLUE) if vals[i] >= 0 else (RED if names[i] not in added else GOLD) for i in range(len(names))]
        ax.barh(range(len(names)), vals, color=cols)
        ax.axvline(0, color=INK, lw=1)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8.5)
        ax.set_xlabel("normalized SVM weight  (mean over CV bins)")
        _caption(ax, title)
        grid(ax, "x")
    added = set(rn) - set(on)
    panel(axes[0][0], rw, "+Prosit (all features)", added)
    if two:
        panel(axes[0][1], ow, "no Prosit (search features only)")
    hs = [Patch(color=BLUE, label="search feature, +weight (→ target)"), Patch(color=RED, label="search feature, −weight"),
          Patch(color=GREEN, label="Prosit-added feature, +weight"), Patch(color=GOLD, label="Prosit-added feature, −weight")]
    fig.legend(handles=hs, loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    _sup(fig, "Percolator feature weights — what drives the score", y=0.995)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    cap = ("Mean normalized Percolator SVM weight per feature (averaged over cross-validation bins), top 22 by "
           "|weight| (exact-zero/unused features dropped). Positive (right) pushes a PSM toward TARGET. Left = the "
           "+Prosit model (green = features Prosit adds: spectral_angle*, pearson_corr, RT, fragment counts…); right = "
           "the no-Prosit model. NB Oktoberfest names the primary search score 'andromeda'/'m0' even for MSFragger.")
    return "weights", "Percolator feature weights", cap, fig_to_uri(fig)


def _qvals(path):
    """Sorted q-values of the target rows in a percolator output file."""
    out = []
    with open(path) as f:
        qi = f.readline().rstrip("\n").split("\t").index("q-value")
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                out.append(float(p[qi]))
            except (IndexError, ValueError):
                pass
    return np.sort(np.asarray(out, float))


def fig_yield(D, P):
    """Identifications accepted as a function of the FDR threshold, +Prosit vs no-Prosit.

    The standard way this comparison is reported in the rescoring literature: one curve per model, so the
    gain is visible at EVERY threshold rather than only at the single 1% operating point.
    """
    lv_have = []
    for lv, name in [("psms", "PSMs"), ("peptides", "Peptides")]:
        r = P["perc"] / f"rescore.percolator.{lv}.txt"
        o = P["perc"] / f"original.percolator.{lv}.txt"
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
        ax.annotate(txt, (FDR, nr), textcoords="offset points", xytext=(12, -6), fontsize=10,
                    color=INK, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#cccccc", alpha=0.9))
        ax.set_xlim(0, 0.05); ax.set_ylim(bottom=0)
        ax.set_xlabel("FDR threshold (q-value)"); ax.set_ylabel(f"accepted target {name} (cumulative)")
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        _caption(ax, f"{name} accepted vs FDR threshold"); grid(ax)
    _sup(fig, "Identification yield vs FDR — with and without Prosit", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    cap = ("Number of target identifications accepted as the FDR threshold is relaxed, for Percolator WITH Prosit "
           "predictions (green) and WITHOUT them (blue); the shaded band between the curves is the gain. The dashed "
           "line marks the 1% operating point, annotated with both counts, the absolute and the relative gain. "
           "Reading: a curve that sits above the other ACROSS the whole range is a genuine discrimination "
           "improvement, not an artefact of where the threshold happens to be drawn.")
    return "yield", "Identification yield vs FDR", cap, fig_to_uri(fig)


def fig_movement(D, P):
    label, score, oscore, q = D["label"], D["score"], D["oscore"], D["q"]
    m = (label == 1) & np.isfinite(score) & np.isfinite(oscore)
    if m.sum() < 50:
        return None
    # accepted sets by PSMId
    resc_acc = ids_at_fdr(P["perc"] / "rescore.percolator.psms.txt", "PSMId")
    orig_acc = ids_at_fdr(P["perc"] / "original.percolator.psms.txt", "PSMId") if (P["perc"] / "original.percolator.psms.txt").exists() else set()
    sid = D["SpecId"]
    in_r = np.array([s in resc_acc for s in sid]); in_o = np.array([s in orig_acc for s in sid])
    gained = m & in_r & ~in_o; lost = m & ~in_r & in_o; common = m & in_r & in_o; neither = m & ~in_r & ~in_o
    fig, ax = plt.subplots(figsize=(8.6, 8))
    rng = np.random.default_rng(0)
    def sub(mask, k=8000):
        idx = np.where(mask)[0]
        return rng.choice(idx, size=min(k, idx.size), replace=False) if idx.size else idx
    ax.scatter(oscore[sub(neither, 6000)], score[sub(neither, 6000)] if False else score[sub(neither, 6000)], s=4, color="#d0d0cc", alpha=0.4, linewidths=0, label=f"neither, n={int(neither.sum()):,}")
    for mask, col, lbl in [(common, COMMON, f"kept (both), n={int(common.sum()):,}"),
                           (gained, GAINED, f"gained by Prosit, n={int(gained.sum()):,}"),
                           (lost, LOST, f"lost by Prosit, n={int(lost.sum()):,}")]:
        s = sub(mask)
        ax.scatter(oscore[s], score[s], s=6, color=col, alpha=0.5, linewidths=0, label=lbl)
    lim = [np.nanpercentile(oscore[m], 0.5), np.nanpercentile(oscore[m], 99.5)]
    ax.plot(lim, lim, color=INK, lw=0.9, ls="--", label="no change")
    ax.set_xlabel("original score (no Prosit)"); ax.set_ylabel("rescored score (+Prosit)")
    ax.legend(frameon=False, fontsize=9, loc="upper left"); _caption(ax, "Per-PSM score movement: original → +Prosit"); grid(ax)
    _sup(fig, "What Prosit moved — original vs rescored score", y=0.99)
    cap = ("Each point is a target PSM: x = its Percolator score without Prosit, y = with Prosit. Colour = acceptance "
           "at 1% FDR — green: accepted only after Prosit (gained), orange: accepted only before (lost), dark-blue: "
           "accepted both, grey: neither. Points above the dashed y=x line were pushed up by Prosit.")
    return "movement", "Rescore-vs-original movement", cap, fig_to_uri(fig)


def fig_per_file(D, P):
    """Identifications PER RAW FILE, +Prosit vs no-Prosit.

    Every other panel in this report pools all raw files together, which hides the question that matters for
    single-cell / low-input runs: does rescoring lift every file, or only rescue a few good ones?
    """
    if "filename" not in D or not np.isfinite(D["oq"]).any():
        return None
    acc_r = (D["label"] == 1) & (D["q"] <= FDR)
    acc_o = (D["label"] == 1) & (D["oq"] <= FDR)
    files = np.array(sorted(set(D["filename"])), dtype=object)
    if files.size < 4:
        return None  # a per-file view needs enough files to be a distribution rather than a list
    pep = np.array([strip_pep(p) for p in D["Peptide"]], dtype=object) if "Peptide" in D else None
    nr, no = [], []
    for f in files:
        inf = D["filename"] == f
        if pep is None:
            nr.append(int((inf & acc_r).sum())); no.append(int((inf & acc_o).sum()))
        else:
            nr.append(len(set(pep[inf & acc_r]))); no.append(len(set(pep[inf & acc_o])))
    nr = np.array(nr, float); no = np.array(no, float)
    unit = "distinct peptides" if pep is not None else "PSMs"
    order = np.argsort(-nr)
    x = np.arange(files.size)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    # (a) the per-file yield curve — dynamic range across files, and the lift on top of it
    ax = axes[0]
    ax.fill_between(x, no[order], nr[order], color=GAINED, alpha=0.18, lw=0, label="gain from Prosit")
    ax.plot(x, no[order], color=COMMON, lw=1.8, label="no Prosit")
    ax.plot(x, nr[order], color=GAINED, lw=1.8, label="+Prosit")
    ax.set_xlabel(f"raw file, ranked by +Prosit yield  (n={files.size} files)")
    ax.set_ylabel(f"{unit} @ 1% FDR"); ax.set_xlim(0, files.size - 1); ax.set_ylim(bottom=0)
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
    n_up = int((nr[ok] > no[ok]).sum()); n_dn = int((nr[ok] < no[ok]).sum())
    from scipy.stats import spearmanr
    r, _ = spearmanr(no[ok], pct)
    ax.text(0.97, 0.94, f"gained in {n_up}/{int(ok.sum())} files · lost in {n_dn}\nSpearman r={r:+.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK,
            bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#cccccc", alpha=0.9))
    ax.set_xlabel(f"{unit} @ 1% FDR without Prosit  (per raw file)")
    ax.set_ylabel("relative gain from Prosit (%)")
    lo = min(0.0, float(pct.min())); hi = float(pct.max())
    ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.18 * (hi - lo))
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    _caption(ax, "(b) relative gain vs how well the file ran"); grid(ax)
    _sup(fig, "Per-raw-file identifications — is the gain uniform?", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    cap = (f"Every other panel pools all raw files; this one splits them. (a) {unit} accepted at 1% FDR per raw file, "
           "files ranked by their +Prosit yield, with vs without Prosit — the spread along x is the run-to-run "
           "variability of the experiment, the shaded band is what rescoring adds on top of it. (b) the same gain "
           "expressed relative to each file's own no-Prosit yield, plotted against that yield: a flat cloud means "
           "rescoring helps every file by roughly the same factor, a downward trend (negative Spearman r) means "
           "weak files benefit disproportionately. NB counts are distinct peptide sequences among PSMs accepted at "
           "1% PSM-level FDR, which is a per-file proxy — Percolator's own peptide-level FDR is estimated globally, "
           "not per file, so these do not sum to the peptide count in the summary.")
    return "perfile", "Per-raw-file identifications", cap, fig_to_uri(fig)


def fig_calibration(D, P):
    from scipy.stats import spearmanr
    label, q = D["label"], D["q"]
    acc = (label == 1) & (q <= FDR); dec = (label == -1)
    have_ard = "abs_rt_diff" in D
    ppm_cols = [(c, lab) for c, lab in [("delta_mass_ppm", "precursor mass error (ppm)"),
                                        ("mean_ppm_error", "fragment mean ppm error (ppm)")] if c in D]
    # NB no pooled observed-vs-predicted iRT scatter here: RT alignment is fitted PER RAW FILE, so pooling
    # all raws bends the band into a curve that says more about the pooling than about the model. The honest
    # per-raw straight-line fits are the iRT-vs-RT SVGs in the per-raw gallery; the pooling-safe summary is
    # the abs_rt_diff residual below.
    panels = []  # list of draw-callables
    if have_ard:
        def _ard(ax):
            va = D["abs_rt_diff"][acc]; vd = D["abs_rt_diff"][dec]
            va = va[np.isfinite(va)]; vd = vd[np.isfinite(vd)]
            hi = np.percentile(np.concatenate([va, vd]), 99); b = np.linspace(0, hi, 60)
            ax.hist(va, bins=b, density=True, histtype="step", lw=2, color=TARGET_C, label=f"accepted targets, n={va.size:,}")
            ax.hist(vd, bins=b, density=True, histtype="step", lw=2, color=DECOY_C, label=f"decoys, n={vd.size:,}")
            ax.axvline(np.median(va), color=TARGET_C, ls=":", lw=1.3); ax.axvline(np.median(vd), color=DECOY_C, ls=":", lw=1.3)
            ax.set_xlabel("|aligned RT − predicted RT|  (abs_rt_diff)"); ax.set_ylabel("density (area=1)")
            ax.legend(frameon=False, fontsize=8.5)
            _caption(ax, f"RT residual (pooling-safe)  median tgt={np.median(va):.2f} vs decoy={np.median(vd):.2f}")
            grid(ax)
        panels.append(_ard)
    for c, lab in ppm_cols:
        def _ppm(ax, c=c, lab=lab):
            va = D[c][acc]; vd = D[c][dec]; va = va[np.isfinite(va)]; vd = vd[np.isfinite(vd)]
            # robust MAD window around the target median: zooms to the monoisotopic core, so isotope-error
            # satellites (±1 Da ≈ hundreds of ppm) and broad decoys spill beyond instead of crushing it.
            med = np.median(va); mad = np.median(np.abs(va - med))
            half = 6 * mad if mad > 1e-6 else max(np.percentile(va, 97) - med, 1.0)
            lo, hi = max(va.min(), med - half), med + half
            if hi - lo < 1e-6:
                hi = lo + 1
            b = np.linspace(lo, hi, 60)
            ax.hist(va, bins=b, range=(lo, hi), density=True, histtype="step", lw=2, color=TARGET_C, label=f"accepted targets, n={va.size:,}")
            ax.hist(vd, bins=b, range=(lo, hi), density=True, histtype="step", lw=2, color=DECOY_C, label=f"decoys, n={vd.size:,}")
            ax.set_xlim(lo, hi); ax.set_xlabel(lab + "  [target core]"); ax.set_ylabel("density (area=1)")
            ax.legend(frameon=False, fontsize=8.5); _caption(ax, lab.split(" (")[0]); grid(ax)
        panels.append(_ppm)
    if not panels:
        return None
    nrow, ncol = ((2, 2) if len(panels) > 2 else (1, len(panels)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 5.2 * nrow), squeeze=False)
    flat = axes.ravel()
    for i, draw in enumerate(panels):
        draw(flat[i])
    for j in range(len(panels), len(flat)):
        flat[j].axis("off")
    _sup(fig, "Calibration — retention time & mass accuracy", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94 if nrow > 1 else 0.9), h_pad=2.5, w_pad=2.5)
    cap = ("Each panel contrasts accepted targets (green) with decoys (grey) — the separation between the two IS the "
           "discriminative power that feature contributes to Percolator. RT panel: the pooling-safe per-PSM aligned RT "
           "residual |aligned RT − predicted RT| (targets tight, decoys broad; medians dotted). The per-raw "
           "observed-vs-predicted iRT fits are in the gallery below, where the per-file alignment is not pooled away. "
           "ppm panels: mass error, zoomed to the TARGET core via a robust MAD window (isotope-error PSMs spill beyond).")
    return "calibration", "RT & mass-error calibration", cap, fig_to_uri(fig)


# ============================== SPECTRA ==============================
GORDER = {"best_score": 0, "cutoff_accepted": 1, "cutoff_rejected": 2, "worst_score": 3,
          "decoy_high": 4, "decoy_low": 5}
TOL_PPM = 20
REPORTER_MIN, REPORTER_MAX = 124.0, 132.0
B_COL, Y_COL, OTH_COL = "#1f6fb2", "#d1352b", "#8a8a85"  # b=blue, y=red, other=grey (mirror-plot convention)
UNM_COL = "#141414"  # unmatched observed peaks: near-black (was light grey) so they read clearly
PROTON = 1.007276


def select_spectra(D, n):
    label, score, q = D["label"], D["score"], D["q"]
    keep = ["SpecId", "filename", "ScanNr", "spectral_angle", "Peptide", "charge",
            "collision_energy_aligned", "RT", "Mass", "score", "q", "label"]
    import pandas as pd
    df = pd.DataFrame({k: D[k] for k in keep if k in D})
    df["pep"] = df.Peptide.map(strip_pep)
    df["ce"] = (df.collision_energy_aligned * 100).round().astype(int)
    df["scan"] = df.ScanNr.astype(int)
    # Mass = theoretical neutral peptide mass (ExpMass in this tab is a bogus sequential index) -> precursor m/z
    df["mz"] = (df.Mass + df.charge * PROTON) / df.charge
    tgt = df[df.label == 1]; dec = df[df.label == -1]
    conf = tgt[tgt.q <= FDR]; rej = tgt[tgt.q > FDR]
    groups = {
        "best_score": tgt.sort_values("score", ascending=False).head(n),
        "worst_score": tgt.sort_values("score", ascending=True).head(n),
        "cutoff_accepted": conf.sort_values("score", ascending=True).head(n),
        "cutoff_rejected": rej.sort_values("score", ascending=False).head(n),
        "decoy_high": dec.sort_values("score", ascending=False).head(n),
        "decoy_low": dec.sort_values("score", ascending=True).head(n),
    }
    rows = []
    for name, g in groups.items():
        gg = g.copy(); gg["group"] = name; rows.append(gg)
    picks = pd.concat(rows, ignore_index=True)
    picks["key"] = list(zip(picks.filename, picks.scan))
    grp_of = picks.groupby("key")["group"].apply(lambda s: sorted(set(s))).to_dict()
    uniq = picks.drop_duplicates("key").reset_index(drop=True)
    return uniq, grp_of


def _ion_type(a):
    return "b" if a and a[0] == "b" else ("y" if a and a[0] == "y" else "o")


def _stems(x, y):
    xs, ys = [], []
    for xi, yi in zip(x, y):
        xs += [xi, xi, None]; ys += [0, yi, None]
    return xs, ys


def build_spectra_fig(D, P, model, n):
    try:
        import pandas as pd
        from pyteomics import mzml as pymzml
        from koinapy import Koina
        import plotly.graph_objects as go
    except Exception as e:
        return None, f"spectra section skipped (missing dep: {e})"
    if P["spectra"] is None:
        return None, "spectra section skipped (no spectra/ dir with mzML)"
    uniq, grp_of = select_spectra(D, n)
    # koina
    kin = pd.DataFrame({"peptide_sequences": uniq.pep.values, "precursor_charges": uniq.charge.values,
                        "collision_energies": uniq.ce.values, "fragmentation_types": ["HCD"] * len(uniq)})
    try:
        pred = Koina(model_name=model, server_url="koina.wilhelmlab.org:443", ssl=True).predict(kin)
    except Exception as e:
        return None, f"spectra section skipped (Koina query failed: {e})"
    _dec = lambda a: a.decode() if isinstance(a, bytes) else str(a)
    lut = {}
    for key, g in pred.groupby(["peptide_sequences", "precursor_charges", "collision_energies"]):
        lut[(str(key[0]), int(key[1]), int(round(float(key[2]))))] = (
            g["mz"].to_numpy(float), g["intensities"].to_numpy(float), [_dec(a) for a in g["annotation"]])

    def pred_row(r):
        return lut.get((str(r.pep), int(r.charge), int(r.ce)), (np.array([]), np.array([]), []))

    # read mzML
    obs_by_key = {}
    for raw, grp in uniq.groupby("filename"):
        mzf = P["spectra"] / f"{raw}.mzML"
        if not mzf.exists():
            continue
        need = {int(s): None for s in grp.scan}
        with pymzml.MzML(str(mzf)) as rd:
            for sp in rd:
                if sp.get("ms level") != 2:
                    continue
                sid = sp.get("id", "")
                sc = int(sid.split("scan=")[-1].split()[0]) if "scan=" in sid else sp.get("index", -1) + 1
                if sc in need:
                    need[sc] = (np.asarray(sp["m/z array"]), np.asarray(sp["intensity array"]))
                if all(v is not None for v in need.values()):
                    break
        for s, v in need.items():
            obs_by_key[(raw, s)] = v

    data = []
    for _, r in uniq.iterrows():
        obs = obs_by_key.get((r.filename, int(r.scan)))
        if obs is None:
            continue
        pmz, pint, pann = pred_row(r)
        data.append(dict(row=r, obs=obs, pmz=pmz, pint=pint, pann=pann,
                         groups=grp_of[(r.filename, int(r.scan))]))
    if not data:
        return None, "spectra section skipped (no observed spectra matched)"

    def gkey(d):
        return (min(GORDER.get(g, 9) for g in d["groups"]), -float(d["row"].score))
    data.sort(key=gkey)

    fig = go.Figure(); spans, annos, xr = [], [], []
    for d in data:
        start = len(fig.data)
        r = d["row"]
        omz = np.array(d["obs"][0]); oin = np.array(d["obs"][1], float)
        rep = (omz >= REPORTER_MIN) & (omz <= REPORTER_MAX); omz, oin = omz[~rep], oin[~rep]
        pmz = np.array(d["pmz"]); pin = np.array(d["pint"], float); pann = d["pann"]
        keep = (pmz > 0) & (pin > 0); pmz, pin, pann = pmz[keep], pin[keep], [a for a, k in zip(pann, keep) if k]
        oin = 100 * oin / oin.max() if oin.size and oin.max() > 0 else oin
        pin = 100 * pin / pin.max() if pin.size and pin.max() > 0 else pin
        otype = np.array(["u"] * len(omz), dtype=object)
        for i, m in enumerate(omz):
            j = np.where(np.abs(pmz - m) <= m * TOL_PPM / 1e6)[0]
            if len(j):
                otype[i] = _ion_type(pann[j[np.argmin(np.abs(pmz[j] - m))]])
        vis0 = len(spans) == 0
        for typ, col, wid in [("u", UNM_COL, 1.0), ("b", B_COL, 1.2), ("y", Y_COL, 1.2), ("o", OTH_COL, 1.2)]:
            sel = otype == typ
            if not sel.any():
                continue
            xs, ys = _stems(omz[sel], oin[sel])
            hov = [f"m/z {m:.4f}<br>{v:.1f}%{' · unmatched' if typ=='u' else ''}" for m, v in zip(omz[sel], oin[sel])]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=col, width=wid),
                                     opacity=0.55 if typ == "u" else 1.0, hoverinfo="text",
                                     text=sum([[h, h, ""] for h in hov], []), visible=vis0, showlegend=False))
        for typ, col in [("b", B_COL), ("y", Y_COL), ("o", OTH_COL)]:
            idx = [i for i, a in enumerate(pann) if _ion_type(a) == typ]
            if not idx:
                continue
            xs, ys = _stems(pmz[idx], -pin[idx])
            hov = [f"{pann[i]}<br>m/z {pmz[i]:.4f}<br>{pin[i]:.1f}%" for i in idx]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=col, width=1.4),
                                     hoverinfo="text", text=sum([[h, h, ""] for h in hov], []),
                                     visible=vis0, showlegend=False))
        spans.append((start, len(fig.data)))
        gs = " / ".join(d["groups"])
        qtxt = "decoy" if int(r.label) == -1 else f"q={float(r.q):.4g}"
        head = (f"<b>{r.pep}</b>  ({gs})<br>"
                f"<span style='font-size:12px'>{r.filename} · scan {int(r.scan)} · z{int(r.charge)} · "
                f"m/z {float(r.mz):.4f} · RT {float(r.RT):.1f} · NCE {int(r.ce)} · "
                f"SA {float(r.spectral_angle):.3f} · score {float(r.score):+.3f} · {qtxt}</span>")
        an = [dict(x=0.0, y=1.10, xref="paper", yref="paper", showarrow=False, xanchor="left",
                   text=head, font=dict(size=13, color=INK))]
        for i in np.argsort(pin)[::-1]:
            if pin.size == 0 or pin[i] < 8:
                break
            an.append(dict(x=float(pmz[i]), y=float(-pin[i] - 4), showarrow=False, text=pann[i],
                           font=dict(size=9, color=B_COL if _ion_type(pann[i]) == "b" else (Y_COL if _ion_type(pann[i]) == "y" else OTH_COL))))
        annos.append(an)
        allmz = np.concatenate([omz, pmz]) if pmz.size else omz
        pad = (allmz.max() - allmz.min()) * 0.03 + 5
        xr.append([float(allmz.min() - pad), float(allmz.max() + pad)])

    buttons = []
    for si, d in enumerate(data):
        vis = [False] * len(fig.data)
        for tt in range(*spans[si]):
            vis[tt] = True
        r = d["row"]; gs = "/".join(d["groups"])
        lab = f"{gs} | {r.pep[:22]} (score {float(r.score):+.2f}, SA {float(r.spectral_angle):.2f})"
        buttons.append(dict(label=lab, method="update", args=[{"visible": vis}, {"annotations": annos[si], "xaxis.range": xr[si]}]))
    for col, nm in [(B_COL, "b ion"), (Y_COL, "y ion"), (UNM_COL, "observed, unmatched")]:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=dict(color=col, width=3), name=nm, showlegend=True))
    fig.update_layout(
        updatemenus=[dict(buttons=buttons, x=0.0, xanchor="left", y=1.20, yanchor="top", showactive=True,
                          direction="down", font=dict(size=11), bgcolor="#f4f4f0", bordercolor="#cccccc")],
        annotations=annos[0], template="plotly_white",
        font=dict(family="Inter, Helvetica, Arial, sans-serif", size=13, color=INK),
        width=1050, height=580, margin=dict(t=150, b=55, l=70, r=40),
        xaxis=dict(title="m/z", range=xr[0], gridcolor="#f2f2ee"),
        yaxis=dict(title="observed ↑ / predicted ↓  (% base peak)", range=[-118, 118], zeroline=False, gridcolor="#eeeeea"),
        legend=dict(orientation="h", y=1.16, x=1.0, xanchor="right", font=dict(size=11)), hovermode="closest")
    fig.add_hline(y=0, line=dict(color="#cccccc", width=0.8))
    html = fig.to_html(full_html=False, include_plotlyjs="inline", div_id="spectra_div")
    note = (f"{len(data)} spectra across groups: " + ", ".join(f"{k}={sum(1 for d in data if k in d['groups'])}" for k in GORDER))
    return html, note


# ============================== native SVGs ==============================
SUMMARY_ORDER = ["psm_1%_FDR", "peptide_1%_FDR", "rescore_original_joint_plot_psm",
                 "rescore_original_joint_plot_peptide", "target_vs_decoys_sa_distribution",
                 "rescore_target_vs_decoys_psm_bins", "original_target_vs_decoys_psm_bins",
                 "rescore_target_vs_decoys_peptide_bins", "original_target_vs_decoys_peptide_bins"]


def collect_svgs(P):
    svgs = {}
    for d in {P["perc"], P["results_dir"]}:
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
    return summary, per_raw


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
    for up in (here.parents[2] if len(here.parents) > 2 else here.parent,
               here.parents[3] if len(here.parents) > 3 else here.parent):
        candidates.append(up / "docs" / "_static" / "img" / "Oktoberfest.svg")
    candidates.append(Path("/cmnfs/home/students/a.weitz/mirror_html_23fp/Oktoberfest.svg"))
    for p in candidates:
        try:
            if p.exists():
                return "data:image/svg+xml;base64," + _b64.b64encode(p.read_bytes()).decode()
        except Exception:
            pass
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
    return (f"<div class='sechead' id='{sid}'><span class='secnum'>{n}</span><h2>{title}</h2></div>"
            f"<div class='secdesc'>{SECTION_DESC.get(sid, '')}</div>")


def _section(n, sid, title, body, is_open=False):
    """A collapsible section: header is the <summary>, content collapses. Summary is open by default."""
    op = " open" if is_open else ""
    desc = SECTION_DESC.get(sid, "")
    return (f"<details class='sec'{op} id='{sid}'><summary>"
            f"<span class='secnum'>{n}</span><span class='sectitle'>{title}</span>"
            f"<span class='secdesc-inline'>{desc}</span></summary>"
            f"<div class='secbody'>{body}</div></details>")


def _kpi(v, l, s="", cls=""):
    sub = f"<div class='s'>{s}</div>" if s else ""
    return f"<div class='kpi {cls}'><div class='v'>{v}</div><div class='l'>{l}</div>{sub}</div>"


def overview_stats(P, D):
    m = (D["label"] == 1) & (D["q"] <= FDR)
    s = dict(n_rows=D["n"], n_t=int((D["label"] == 1).sum()), n_d=int((D["label"] == -1).sum()),
             median_sa=float(np.median(D["spectral_angle"][m])) if m.any() else float("nan"),
             cutoff=float(D["score"][m].min()) if m.any() else float("nan"))
    rp = P["perc"]
    rid = ids_at_fdr(rp / "rescore.percolator.psms.txt", "PSMId"); s["acc_psm"] = len(rid)
    pep_f = rp / "rescore.percolator.peptides.txt"
    rpep = ids_at_fdr(pep_f, "peptide") if pep_f.exists() else set(); s["acc_pep"] = len(rpep)
    op, opep = rp / "original.percolator.psms.txt", rp / "original.percolator.peptides.txt"
    if op.exists():
        oid = ids_at_fdr(op, "PSMId")
        s.update(gain_psm=len(rid) - len(oid), gained_psm=len(rid - oid), lost_psm=len(oid - rid))
    if opep.exists() and rpep:
        oidp = ids_at_fdr(opep, "peptide")
        s.update(gain_pep=len(rpep) - len(oidp), gained_pep=len(rpep - oidp), lost_pep=len(oidp - rpep))
    return s


def kpi_cards(s):
    c = [_kpi(f"{s['acc_psm']:,}", "Target PSMs @1% FDR", f"of {s['n_t']:,} target PSMs"),
         _kpi(f"{s['acc_pep']:,}", "Peptides @1% FDR", "")]
    if "gain_psm" in s:
        c.append(_kpi(f"{s['gain_psm']:+,}", "Prosit net PSMs", f"+{s['gained_psm']:,} / −{s['lost_psm']:,}",
                      "pos" if s["gain_psm"] >= 0 else "neg"))
    if "gain_pep" in s:
        c.append(_kpi(f"{s['gain_pep']:+,}", "Prosit net peptides", f"+{s['gained_pep']:,} / −{s['lost_pep']:,}",
                      "pos" if s["gain_pep"] >= 0 else "neg"))
    c += [_kpi(f"{s['median_sa']:.3f}", "Median SA (accepted)"),
          _kpi(f"{s['n_d']:,}", "Decoy PSMs"),
          _kpi(f"{s['cutoff']:.2f}", "Score at 1% cutoff")]
    return "<div class='kpis'>" + "".join(c) + "</div>"


def config_rows(cfg):
    g = cfg.get; models = cfg.get("models", {}) or {}; inp = cfg.get("inputs", {}) or {}
    rows = [("Workflow", g("type")), ("Search engine", inp.get("search_results_type")),
            ("Spectra type", inp.get("spectra_type")), ("Intensity model", models.get("intensity")),
            ("iRT model", models.get("irt")), ("Prediction server", g("prediction_server")),
            ("Mass tolerance", f"{g('massTolerance')} {g('unitMassTolerance', '')}".strip() if g("massTolerance") is not None else None),
            ("FDR estimation", g("fdr_estimation_method")), ("Peak matching", g("matching_method")),
            ("All features", str(g("allFeatures")) if g("allFeatures") is not None else None),
            ("TMT tag", g("tag"))]
    return [(k, v) for k, v in rows if v not in (None, "", "None", "None ")]


def build_report_html(P, D, model, run_name, cfg, stats, sections, spectra_html, spectra_note,
                      summary_svgs, per_raw, want_spectra):
    nav = [("summary", "Summary")] + [(sid, title) for sid, title, _, _ in sections]
    if spectra_html:
        nav.append(("spectra", "Spectra"))
    if summary_svgs:
        nav.append(("native", "Native plots"))
    if per_raw:
        nav.append(("gallery", "Per-raw gallery"))
    engine = (cfg.get("inputs", {}) or {}).get("search_results_type", "")
    badges = "".join(f"<span class='badge'>{b}</span>" for b in [model] + ([engine] if engine else []))
    badges += f"<span class='badge grey'>generated {datetime.now():%Y-%m-%d %H:%M}</span>"

    h = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width,initial-scale=1'>",
         f"<title>Oktoberfest report — {run_name}</title><style>{CSS}</style></head><body><div class='wrap'>"]
    logo = _logo_uri()
    logo_html = f"<img class='logo' src='{logo}' alt='Oktoberfest'>" if logo else ""
    h.append(f"<div class='masthead'>{logo_html}<div class='kicker'>Rescoring report</div>"
             f"<h1>{run_name}</h1><div class='run'>Spectral-prediction PSM rescoring — diagnostic report</div>"
             f"<div class='badges'>{badges}</div></div>")
    h.append("<nav class='toc'>" + "".join(f"<a href='#{sid}'>{t}</a>" for sid, t in nav) + "</nav>")

    n = 1
    rows = config_rows(cfg)
    cfg_html = ("<div class='cfg'><table>" + "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows) + "</table></div>") if rows else ""
    h.append(_section(n, "summary", "Summary", kpi_cards(stats) + cfg_html, is_open=True))  # only Summary open
    for sid, title, cap, uri in sections:
        n += 1
        body = (f"<div class='figcard'><img src='{uri}' alt='{title}'>"
                f"<div class='figcap'><b class='h'>What this shows</b>{cap}</div></div>")
        h.append(_section(n, sid, title, body))
    if spectra_html:
        n += 1
        legend = ("<div class='legendbar'>"
                  f"<span><span class='sw' style='background:{B_COL}'></span><b>b</b> ion (matched)</span>"
                  f"<span><span class='sw' style='background:{Y_COL}'></span><b>y</b> ion (matched)</span>"
                  "<span><span class='sw' style='background:#141414'></span>observed, unmatched</span>"
                  "<span>observed ↑ / predicted ↓ · TMT reporter 124–132 m/z removed · 20 ppm match window</span></div>")
        body = ("<div class='spectra-wrap'>" + legend +
                "<div class='figcap' style='margin:0 0 8px'><b class='h'>How to read</b>"
                "Observed spectrum (up, from mzML) mirrored against the clean Koina re-query (down). "
                "Use the dropdown to switch PSM. Each header carries the full info line: "
                "raw · scan · peptide(+mods) · charge · m/z · RT · NCE · SA · score · q. Groups: best/worst percolator "
                "score, both sides of the score≈0 cutoff, plus high- and low-scoring decoys.</div>"
                f"<div class='small' style='margin:0 0 8px'>{spectra_note}</div>" + spectra_html + "</div>")
        h.append(_section(n, "spectra", "Spectra viewer", body))
    elif want_spectra:
        n += 1
        h.append(_section(n, "spectra", "Spectra viewer", f"<div class='figcard'><div class='warn'>{spectra_note}</div></div>"))
    if summary_svgs:
        n += 1
        half = {"psm_1%_FDR", "peptide_1%_FDR"}  # count plots -> render at half size
        cells = "".join(f"<div class='svgcard{' half' if stem in half else ''}'><div class='t'>{stem}</div>"
                        f"<img src='{svg_to_uri(p)}'></div>" for stem, p in summary_svgs)
        h.append(_section(n, "native", "Native Oktoberfest plots", f"<div class='grid2'>{cells}</div>"))
    if per_raw:
        n += 1
        n_full = sum(1 for r in per_raw.values() if len(r) >= 3)
        note = ""
        if n_full < len(per_raw):
            note = (f"<div class='secdesc' style='margin:0 0 10px 0'>Note: Oktoberfest emitted the iRT-vs-RT plot for "
                    f"only {n_full}/{len(per_raw)} raw files, so some entries show 2 plots (CE violin + mean) rather "
                    f"than 3 — the missing SVG was never produced by the run.</div>")
        body = note + "".join(
            f"<details class='raw'><summary>{raw} · {len(per_raw[raw])} plots</summary>"
            f"<div class='grid2' style='margin:8px 0'>" +
            "".join(f"<div class='svgcard'><div class='t'>{k}</div><img src='{svg_to_uri(v)}'></div>"
                    for k, v in sorted(per_raw[raw].items())) + "</div></details>"
            for raw in sorted(per_raw))
        h.append(_section(n, "gallery", "Per-raw gallery", body))
    h.append(f"<footer>Generated by <b>oktoberfest_investigate.py</b> · source: {P['perc']} · "
             f"{datetime.now():%Y-%m-%d %H:%M}.<br>Non-native panels define every colour, mark, unit and n in-figure. "
             "Spectra are observed mzML vs a clean Koina re-query (not the rescoring's stored predictions). "
             "A diagnostic aid alongside — not a replacement for — the primary Oktoberfest outputs.</footer>")
    # open a collapsed section when its nav link (or any in-page anchor) targets it
    h.append("<script>function _openHash(){var el=document.getElementById((location.hash||'').slice(1));"
             "if(el&&el.tagName==='DETAILS'){el.open=true;el.scrollIntoView();}}"
             "window.addEventListener('hashchange',_openHash);window.addEventListener('load',_openHash);</script>")
    h.append("</div></body></html>")
    return "\n".join(h)


def build_report(out_dir, out_html=None, n_per_group=20, want_spectra=True, log=print):
    """Build the HTML report for an Oktoberfest run. `out_dir` may be the run's out/ folder or a
    percolator dir. Returns the written path. Individual figures/spectra are already guarded; callers
    that must never fail (e.g. inside the pipeline) should still wrap this in try/except (see
    build_report_safe)."""
    P = find_paths(out_dir)
    cfg = {}
    model = "Prosit_2025_intensity_40PTM"
    if P["config"]:
        try:
            cfg = json.loads(Path(P["config"]).read_text())
            model = cfg.get("models", {}).get("intensity", model)
        except Exception:
            cfg = {}
    run_name = P["run_dir"].name
    engine = (cfg.get("inputs", {}) or {}).get("search_results_type", "")
    RUN_CTX["suffix"] = "  ·  ".join(x for x in [engine, model, run_name] if x)  # adaptive figure subtitles
    out_html = Path(out_html) if out_html else (P["run_dir"] / f"investigate_{run_name}.html")

    log(f"[investigate] run={run_name} perc={P['perc']} model={model}")
    D = load_data(P)

    sections = []
    for fn in [fig_yield, fig_movement, fig_per_file, fig_weights, fig_sa_features, fig_headroom, fig_calibration]:
        try:
            res = fn(D, P)
        except Exception as e:
            log(f"[investigate]   skip {fn.__name__}: {e}"); res = None
        if res:
            sections.append(res)

    spectra_html, spectra_note = None, "spectra section disabled"
    if want_spectra:
        try:
            spectra_html, spectra_note = build_spectra_fig(D, P, model, n_per_group)
        except Exception as e:
            spectra_html, spectra_note = None, f"spectra section skipped ({e})"
        log(f"[investigate]   {spectra_note}")

    summary_svgs, per_raw = collect_svgs(P)
    stats = overview_stats(P, D)
    html = build_report_html(P, D, model, run_name, cfg, stats, sections, spectra_html, spectra_note,
                             summary_svgs, per_raw, want_spectra)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html)
    log(f"[investigate] wrote {out_html} ({out_html.stat().st_size / 1e6:.1f} MB, {len(sections)} figs, "
        f"{'spectra' if spectra_html else 'no spectra'}, {len(per_raw)} raw galleries)")
    return out_html


def build_report_safe(out_dir, out_html=None, n_per_group=20, want_spectra=True, log=print):
    """Never-raises wrapper for pipeline use: catches EVERYTHING and returns the path or None."""
    try:
        return build_report(out_dir, out_html=out_html, n_per_group=n_per_group,
                             want_spectra=want_spectra, log=log)
    except Exception as e:  # noqa: BLE001 - deliberately swallow all errors so a report never breaks a run
        try:
            log(f"[investigate] report generation failed, skipped: {e}")
        except Exception:
            pass
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("out_html", nargs="?", default=None)
    ap.add_argument("--n-per-group", type=int, default=20)
    ap.add_argument("--no-spectra", action="store_true")
    a = ap.parse_args()
    build_report(a.out_dir, a.out_html, a.n_per_group, want_spectra=not a.no_spectra)


if __name__ == "__main__":
    main()
