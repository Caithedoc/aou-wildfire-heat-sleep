#!/usr/bin/env python3
"""JNO Supplement 1 eFigures (2026-09-22), same visual system as the main figures (make_fig_jno_v9.py):
no in-plot titles, Helvetica, true minus signs, italic P without leading zero.
Numbers are read from the eTables of the current Supplement docx (single source of truth), except
eFigure 6, which uses the aggregate CSVs exported from the Workbench (rev_fig1_timeline/region.csv).

eFigure 1  Participant and night flowchart, black and white (eTable 1)  [表号按 2026-09-24 重编号: dlag=eTable 3, sensitivity=5, subgroups=6, causal=7, neg controls=9]
eFigure 2  Subgroup associations, marginal (A) and joint-exposure (B) (eTable 6)
eFigure 3  Distributed lag, lags 0-7 (eTable 4)
eFigure 4  Causal-forest calibration (eTable 7)
eFigure 5  Marginal smoke-night association across specifications (eTables 3 and 8)
eFigure 6  Exposure nights over time (A) and by region (B) (was panels B-C of the old overview figure)
"""
import os
import re


def rng(s):
    """Only numeric ranges get an en dash; units get proper superscripts."""
    s = re.sub(r"(\d)-(\d)", "\\1\u2013\\2", s)
    s = re.sub(r"(\dxx)-(\dxx)", "\\1\u2013\\2", s)
    return s.replace("µg/m3", "µg/m³")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.oxml.ns import qn
from matplotlib.patches import FancyBboxPatch
from lu_palette import CAT as LU_CAT, ACC as LU_ACC, LU_BLUE, LU_RED, BAND_ALPHA

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)
SUPP = f"{BASE}/submission_JNO/CZ_SmokeHeat_Sleep_Supplement1_JNO_20260922.docx"
OUT = f"{BASE}/figures"
INK, INK2, GRID = "#1a1a1a", "#555555", "#d9d9d9"
ACC = LU_ACC  # Lu 2025 red accent
CAT = LU_CAT  # Lu 2025 categorical: orange / blue / red

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 8.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.7, "axes.edgecolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
    "axes.unicode_minus": True, "svg.fonttype": "none", "pdf.fonttype": 42,
})


def m(x, nd=2):
    return f"{x:.{nd}f}".replace("-", "−")


def num(s):
    s = s.strip().replace("−", "-").replace(",", "")
    return float(s)


def ci(s):
    a, b = s.split(" to ")
    return num(a), num(b)


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=600, bbox_inches="tight",
                    metadata={"Creator": "Zhang, Cai"} if ext != "png" else None)
    plt.close(fig)
    print("saved", name)


# ---------------------------------------------------------------- read eTables from the docx
doc = Document(SUPP)
body = list(doc.element.body)
tables = {}
for i, el in enumerate(body[:-1]):
    t = "".join(x.text or "" for x in el.iter(qn("w:t")))
    mm = re.match(r"eTable (\d+)\.", t)
    if mm and el.find(".//" + qn("w:tab")) is None and body[i + 1].tag == qn("w:tbl"):
        rows = []
        for tr in body[i + 1].iter(qn("w:tr")):
            rows.append(["".join(x.text or "" for x in tc.iter(qn("w:t"))) for tc in tr.findall(qn("w:tc"))])
        tables[int(mm.group(1))] = rows

# ================================================================ eFigure 1: flowchart (B/W)
t1 = tables[1][1:]
steps = [(r[0].split(". ", 1)[1], r[1], r[2]) for r in t1]
excl = [None] + [(r[3], r[4]) for r in t1[1:]]
MAIN_LABELS = [
    "Participants with main sleep records\n(Controlled Tier CDR version 9)",
    "Nights from 2016 to 2023",
    "Nights passing cleaning criteria with a\nvalid 3-digit ZIP code",
    "Person-months linked to daily exposures\ncontaining at least 1 smoke night",
    "Analytic sample: strata with at least 2 nights\nand complete covariates",
]
EXCL_LABELS = [
    None,
    "Outside the 2016-2023 exposure window",
    "Failed nightly cleaning (duplicate records;\ntotal sleep time outside 180-840 min)\nor no valid 3-digit ZIP code",
    "Person-months without a smoke night\n(uninformative for the self-matched design)",
    "Strata with fewer than 2 nights or\nincomplete covariates",
]
fig, ax = plt.subplots(figsize=(7.0, 7.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
H, W, XC = 1.12, 4.6, 0.3
ys = [9.0 - k * 2.0 for k in range(5)]


def box(x, y, w, h, text, lw=0.9, bold_last=True):
    ax.add_patch(FancyBboxPatch((x, y - h / 2), w, h, boxstyle="square,pad=0", fc="white", ec=INK, lw=lw))
    ax.text(x + 0.18, y, text, ha="left", va="center", fontsize=7.8, color=INK, linespacing=1.35)


for k, (y, lab, (name, npart, nnight)) in enumerate(zip(ys, MAIN_LABELS, steps)):
    txt = rng(f"{lab}\n{npart} participants; {nnight} nights")
    box(XC, y, W, H, txt, lw=1.4 if k == len(ys) - 1 else 0.9)
    if k < len(ys) - 1:
        ax.annotate("", xy=(XC + W / 2, ys[k + 1] + H / 2), xytext=(XC + W / 2, y - H / 2),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9, mutation_scale=9))
for k in range(1, 5):
    yk = (ys[k - 1] + ys[k]) / 2
    ex_p, ex_n = excl[k]
    ax.plot([XC + W / 2, 5.55], [yk, yk], color=INK, lw=0.9)
    ax.annotate("", xy=(5.55, yk), xytext=(5.45, yk),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9, mutation_scale=9))
    txt = rng(f"Excluded: {ex_p} participants; {ex_n} nights\n{EXCL_LABELS[k]}")
    ax.add_patch(FancyBboxPatch((5.55, yk - 0.62), 4.35, 1.24, boxstyle="square,pad=0", fc="white",
                                ec=INK, lw=0.7, ls=(0, (3, 2))))
    ax.text(5.7, yk, txt, ha="left", va="center", fontsize=7.2, color=INK, linespacing=1.3)
save(fig, "eFig_JNO1_flowchart")

# ================================================================ eFigure 2: subgroups (A marginal, B joint)
t6 = tables[6]
hdr = t6[0]
rowsA, rowsB, cur = [], [], None
for r in t6[1:]:
    if len(r) == 1 or r[1] == r[0]:
        cur = rowsA if "Panel A" in r[0] else rowsB
        continue
    b = num(r[1]); lo, hi = ci(r[2])
    pint = r[-1].strip()
    cur.append((r[0], b, lo, hi, pint))
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(wspace=0.95))
for ax, rows, lab, xl, pc in ((axes[0], rowsA, "A", (-3.0, 1.2), LU_BLUE), (axes[1], rowsB, "B", (-6.0, 1.5), LU_RED)):
    y = np.arange(len(rows))[::-1]
    ax.axvline(0, color=INK2, lw=0.7, ls="--", zorder=1)
    for yy, (name, b, lo, hi, pint) in zip(y, rows):
        ax.plot([lo, hi], [yy, yy], color=pc, lw=1.5, solid_capstyle="round", zorder=2, alpha=0.85)
        ax.scatter([b], [yy], s=28, color=pc, zorder=3, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7.8)
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlim(*xl)
    ax.set_xlabel("Beta, min" if lab == "A" else "Joint-exposure night vs neither, min")
    ax.text(-0.02, 1.03, lab, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")
axes[0].set_xlabel("Smoke night vs non-smoke night, min")
save(fig, "eFig_JNO2_subgroups")

# ================================================================ eFigure 3: distributed lag
t4 = tables[3][1:]  # eTable 3 (dlag) after 2026-09-24 renumbering
lag = np.array([int(r[0]) for r in t4]); b = np.array([num(r[1]) for r in t4])
lo = np.array([ci(r[2])[0] for r in t4]); hi = np.array([ci(r[2])[1] for r in t4])
fig, ax = plt.subplots(figsize=(5.6, 3.0))
ax.axhline(0, color=INK2, lw=0.7, ls="--")
ax.fill_between(lag, lo, hi, color=LU_RED, alpha=BAND_ALPHA, lw=0)
ax.plot(lag, b, "-o", color=LU_RED, lw=1.5, ms=5, mec="white", mew=0.8)
ax.set_xticks(lag)
ax.set_xlabel("Days after smoke night (lag)")
ax.set_ylabel("Change in total sleep time, min")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
save(fig, "eFig_JNO3_distributed_lag")

# ================================================================ eFigure 4: causal-forest calibration
t7 = tables[7][1:]
q = np.arange(1, len(t7) + 1)
pred = [num(r[1]) for r in t7]; obs = [num(r[2]) for r in t7]
olo = [ci(r[3])[0] for r in t7]; ohi = [ci(r[3])[1] for r in t7]
fig, ax = plt.subplots(figsize=(5.0, 3.3))
ax.axhline(0, color=INK2, lw=0.7, ls="--")
ax.plot(q, pred, "-s", color=LU_BLUE, lw=1.1, ms=5, mfc="white", mec=LU_BLUE, label="Predicted by causal forest")
for i in range(len(q)):
    ax.plot([q[i], q[i]], [olo[i], ohi[i]], color=ACC, lw=1.3)
ax.plot(q, obs, "-o", color=ACC, lw=1.5, ms=5, mec="white", mew=0.8, label="Observed (sorted-group effect)")
ax.set_xticks(q); ax.set_xticklabels(["1\nMost\nnegative", "2", "3", "4", "5\nMost\npositive"])
ax.set_xlabel("Quintile of predicted individual effect")
ax.set_ylabel("Change in total sleep time, min")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
save(fig, "eFig_JNO4_causal_forest")

# ================================================================ eFigure 5: robustness forest
t3 = tables[5][1:]; t8 = tables[10][1:]  # eTable 5 sensitivity, eTable 10 negative controls
spec = [(r[0], num(r[1]), *ci(r[2]), "spec") for r in t3]
spec += [("Negative control: " + r[0][0].lower() + r[0][1:], num(r[1]), *ci(r[2]), "neg") for r in t8]
fig, ax = plt.subplots(figsize=(7.0, 5.2))
y = np.arange(len(spec))[::-1]
ax.axvline(0, color=INK2, lw=0.7, ls="--", zorder=1)
for yy, (name, b, lo, hi, kind) in zip(y, spec):
    primary = name.startswith("Primary")
    col = ACC if primary else INK
    ax.plot([lo, hi], [yy, yy], color=col if kind != "neg" else INK2, lw=1.4, solid_capstyle="round", zorder=2)
    if kind == "neg":
        ax.scatter([b], [yy], s=30, facecolor="white", edgecolor=INK, linewidth=1.0, zorder=3)
    else:
        ax.scatter([b], [yy], s=30, color=col, edgecolor="white", linewidth=0.8, zorder=3)
    ax.text(3.05, yy, f"{m(b)} ({m(lo)} to {m(hi)})", va="center", ha="left", fontsize=7.2, color=INK)
ax.set_yticks(y); ax.set_yticklabels([rng(s[0]) for s in spec], fontsize=7.6)
ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
ax.set_xlim(-2.0, 3.0); ax.spines["bottom"].set_bounds(-2, 3)
ax.set_xlabel("Change in total sleep time on smoke nights, min")
ax.text(3.05, len(spec) - 0.3, "Beta, min (95% CI)", fontsize=7.4, fontweight="bold", color=INK, va="bottom")
save(fig, "eFig_JNO5_robustness")

# ================================================================ eFigure 6: exposure nights over time and by region
tl = pd.read_csv(f"{BASE}/results_import/rev_fig1_timeline.csv")
rg = pd.read_csv(f"{BASE}/results_import/rev_fig1_region.csv")
KINDS = [("hot_only", "Hot only"), ("smoke_only", "Smoky only"), ("both", "Hot and smoky")]
fig = plt.figure(figsize=(7.2, 5.6))
axA = fig.add_axes([0.1, 0.52, 0.86, 0.42])
tlp = tl.dropna(subset=["n"]).copy()
tlp["t"] = tlp["year"] + (tlp["month"] - 0.5) / 12
piv = (tlp[tlp.kind != "neither"].pivot_table(index="t", columns="kind", values="n", aggfunc="sum")
       .reindex(columns=[k for k, _ in KINDS]).fillna(0).sort_index())
bottom = np.zeros(len(piv))
for k, lab in KINDS:
    axA.bar(piv.index, piv[k] / 1000, width=1 / 13, bottom=bottom / 1000, color=CAT[k], label=lab,
            lw=0.3, edgecolor="white")
    bottom += piv[k].values
axA.set_ylabel("Exposure nights per month, thousands")
axA.set_xlim(2015.9, 2024.1); axA.set_xticks(range(2016, 2024))
axA.yaxis.grid(True, color=GRID, lw=0.5); axA.set_axisbelow(True)
axA.legend(frameon=False, fontsize=7.5, loc="upper left")
axA.text(-0.02, 1.03, "A", transform=axA.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")

axB = fig.add_axes([0.25, 0.07, 0.71, 0.3])
rgp = rg.dropna(subset=["n"]).copy()
rgp["region"] = np.select([rgp.zip1.astype(str).isin(list("012")), rgp.zip1.astype(str).isin(list("34567")),
                           rgp.zip1.astype(str).isin(list("89"))],
                          ["East (ZIP 0-2)", "Central (ZIP 3-7)", "West (ZIP 8-9)"], default="other")
order = ["East (ZIP 0-2)", "Central (ZIP 3-7)", "West (ZIP 8-9)"]
agg = (rgp[rgp.region != "other"].pivot_table(index="region", columns="kind", values="n", aggfunc="sum")
       .reindex(index=order).reindex(columns=[k for k, _ in KINDS]).fillna(0))
yp = np.arange(len(agg))[::-1]; left = np.zeros(len(agg))
for k, lab in KINDS:
    axB.barh(yp, agg[k] / 1000, left=left / 1000, color=CAT[k], height=0.55, lw=0.3, edgecolor="white")
    left += agg[k].values
for yy, (_, row) in zip(yp, agg.iterrows()):
    axB.text(row.sum() / 1000 + 2, yy, f"{int(row['both']):,} hot and smoky", fontsize=7.3, va="center", color=INK)
axB.set_yticks(yp); axB.set_yticklabels([o.replace("-", "–") for o in order])
axB.tick_params(axis="y", length=0); axB.spines["left"].set_visible(False)
axB.set_xlabel("Exposure nights, thousands")
axB.set_xlim(0, left.max() / 1000 * 1.28)
axB.text(-0.30, 1.05, "B", transform=axB.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")
save(fig, "eFig_JNO6_exposure_time_region")
