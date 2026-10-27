#!/usr/bin/env python3
"""New Figure 1: study overview. A: participant/night flow (attrition S1-S6);
B: monthly exposure-night counts 2016-2023 by type; C: exposure mix by ZIP region.
Data: attrition (published Table S1 numbers) + results_import/rev_fig1_timeline.csv
+ rev_fig1_region.csv. Outputs Fig1_overview.png/svg/pdf."""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)
tl = pd.read_csv(f"{BASE}/results_import/rev_fig1_timeline.csv")
rg = pd.read_csv(f"{BASE}/results_import/rev_fig1_region.csv")

plt.rcParams.update({"font.family": "Helvetica", "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8})

C_HOT, C_SMOKE, C_BOTH = "#e8a838", "#7f7f7f", "#b2182b"

fig = plt.figure(figsize=(7.2, 8.6))

# ---------------------------------------------------------------- A: flow
axA = fig.add_axes([0.06, 0.66, 0.88, 0.30]); axA.axis("off")
steps = [
    ("All of Us participants with main-sleep records (CDR v9)", "62,911 participants · 39,303,369 nights"),
    ("Restricted to study window 2016–2023", "58,885 participants · 31,813,502 nights"),
    ("Night-level cleaning (deduplication; 3–14 h; efficiency ≤100%)", "58,262 participants · 29,957,387 nights"),
    ("Linked to 3-digit ZIP code and daily exposures;\nperson-months containing a smoke night", "31,416 participants · 1,669,021 nights"),
    ("Strata with ≥2 eligible nights and complete covariates\n(regression sample)", "30,268 participants · 1,633,414 nights"),
]
n = len(steps)
for i, (lab, num) in enumerate(steps):
    y1 = 1 - i * (1.0 / n)
    axA.add_patch(plt.Rectangle((0.03, y1 - 0.155), 0.94, 0.15, fill=True,
                                facecolor="#f5f0e8" if i < n - 1 else "#fbe3c9",
                                edgecolor="0.55", lw=0.7))
    axA.text(0.055, y1 - 0.08, lab, fontsize=7.6, va="center")
    axA.text(0.965, y1 - 0.08, num, fontsize=7.6, va="center", ha="right",
             fontweight="bold")
    if i < n - 1:
        axA.annotate("", xy=(0.5, y1 - 0.195), xytext=(0.5, y1 - 0.157),
                     arrowprops=dict(arrowstyle="-|>", color="0.4", lw=0.9))
axA.set_xlim(0, 1); axA.set_ylim(-0.05, 1.02)
axA.set_title("A   Cohort construction", loc="left", fontsize=10, fontweight="bold")

# ---------------------------------------------------------------- B: timeline
axB = fig.add_axes([0.10, 0.375, 0.86, 0.225])
tlp = tl.dropna(subset=["n"]).copy()
tlp["t"] = tlp["year"] + (tlp["month"] - 0.5) / 12
piv = (tlp[tlp.kind != "neither"]
       .pivot_table(index="t", columns="kind", values="n", aggfunc="sum")
       .reindex(columns=["hot_only", "smoke_only", "both"]).fillna(0).sort_index())
bottom = np.zeros(len(piv))
w = 1 / 13
for kind, color, lab in [("hot_only", C_HOT, "Hot only"),
                         ("smoke_only", C_SMOKE, "Smoky only"),
                         ("both", C_BOTH, "Hot + smoky")]:
    axB.bar(piv.index, piv[kind] / 1000, width=w, bottom=bottom / 1000,
            color=color, label=lab, lw=0)
    bottom += piv[kind].values
axB.set_ylabel("Exposure nights per month\n(thousands)")
axB.set_xlim(2015.9, 2024.1)
axB.set_xticks(range(2016, 2024))
axB.legend(frameon=False, fontsize=7.5, loc="upper left")
axB.set_title("B   Exposure nights over time (2023 contributes 57% of analytic nights)",
              loc="left", fontsize=10, fontweight="bold")

# ---------------------------------------------------------------- C: region
axC = fig.add_axes([0.10, 0.045, 0.86, 0.24])
rgp = rg.dropna(subset=["n"]).copy()
rgp["region"] = np.select(
    [rgp.zip1.astype(str).isin(list("012")),
     rgp.zip1.astype(str).isin(list("34567")),
     rgp.zip1.astype(str).isin(list("89"))],
    ["East (ZIP 0xx–2xx)", "Central/South (ZIP 3xx–7xx)", "West (ZIP 8xx–9xx)"],
    default="other")
agg = (rgp[rgp.region != "other"]
       .pivot_table(index="region", columns="kind", values="n", aggfunc="sum")
       .reindex(index=["East (ZIP 0xx–2xx)", "Central/South (ZIP 3xx–7xx)",
                       "West (ZIP 8xx–9xx)"])
       .reindex(columns=["hot_only", "smoke_only", "both"]).fillna(0))
ypos = np.arange(len(agg))[::-1]
left = np.zeros(len(agg))
for kind, color, lab in [("hot_only", C_HOT, "Hot only"),
                         ("smoke_only", C_SMOKE, "Smoky only"),
                         ("both", C_BOTH, "Hot + smoky")]:
    axC.barh(ypos, agg[kind] / 1000, left=left / 1000, color=color, height=0.55, lw=0)
    left += agg[kind].values
for y, (idx, row) in zip(ypos, agg.iterrows()):
    tot = row.sum()
    axC.text(tot / 1000 + 3, y, f"both: {int(row['both']):,}", fontsize=7.5,
             va="center", color=C_BOTH)
axC.set_yticks(ypos)
axC.set_yticklabels(agg.index, fontsize=8.5)
axC.tick_params(axis="y", length=0)
axC.set_xlabel("Exposure nights (thousands)")
axC.set_title("C   Exposure nights by ZIP-code region", loc="left",
              fontsize=10, fontweight="bold")
axC.spines["left"].set_visible(False)

for ext, dpi in (("png", 600), ("svg", None), ("pdf", None)):
    fig.savefig(f"{BASE}/figures/Fig1_overview.{ext}",
                dpi=dpi if dpi else "figure", bbox_inches="tight")
print("saved Fig1_overview.*")
