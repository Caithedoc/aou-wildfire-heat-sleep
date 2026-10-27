#!/usr/bin/env python3
"""New Figure 2: continuous joint smoke x temperature response curves.
Data: results_import/rev_joint_curves.csv (VM spline FE model, delta-method CIs).
Outputs PNG (600 dpi) + SVG + PDF, new numbering Fig2_joint_response_curves.*"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)
df = pd.read_csv(f"{BASE}/results_import/rev_joint_curves.csv")

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
})

SPEC = [
    (0.50, "#3b7bbf", "50th (median night)"),
    (0.75, "#e8a838", "75th"),
    (0.90, "#e2662e", "90th (hot-night threshold)"),
    (0.95, "#b2182b", "95th"),
]

fig, ax = plt.subplots(figsize=(7.0, 4.6))
for tp, color, lab in SPEC:
    d = df[df.tpct == tp]
    ax.plot(d.smoke, d.dTST, color=color, lw=2.0, label=lab, zorder=3)
    ax.fill_between(d.smoke, d.lo, d.hi, color=color, alpha=0.10, lw=0, zorder=1)

ax.axhline(0, color="0.4", lw=0.8, ls="--", zorder=0)
for thr in (20,):
    ax.axvline(thr, color="0.75", lw=0.7, ls=":", zorder=0)
ax.text(20, ax.get_ylim()[0], "", fontsize=7)

ax.set_xlabel("Wildfire-specific PM$_{2.5}$ (µg/m$^3$)")
ax.set_ylabel("Change in total sleep time vs. no smoke (minutes)")
ax.set_title("Smoke-related sleep loss grows with same-night temperature",
             fontsize=10.5, pad=10)
ax.set_xlim(0, 60)
ax.set_ylim(-3.4, 0.4)

leg = ax.legend(title="Temperature percentile\n(within ZIP3–calendar month)",
                frameon=False, fontsize=8, title_fontsize=8, loc="lower left")
ax.annotate("smoke × temperature interaction p = 0.011\n(spline joint model)",
            xy=(0.985, 0.955), xycoords="axes fraction",
            ha="right", va="top", fontsize=7.5, color="0.35")
ax.annotate("primary smoke-night\nthreshold (20 µg/m$^3$)",
            xy=(20, 0.28), ha="center", va="bottom", fontsize=7, color="0.55")

fig.tight_layout()
for ext, dpi in (("png", 600), ("svg", None), ("pdf", None)):
    fig.savefig(f"{BASE}/figures/Fig2_joint_response_curves.{ext}",
                dpi=dpi if dpi else "figure", bbox_inches="tight")
print("saved Fig2_joint_response_curves.png/svg/pdf")
