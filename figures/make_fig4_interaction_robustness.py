#!/usr/bin/env python3
"""New Figure 4: robustness of the smoke-heat interaction.
Panel A: temperature-specification ladder x variance-clustering level
         (betas/CIs from rev2_fourcell_ladder excess rows + interaction ladder;
         values transcribed from results_import/rev2_key_extracts.txt).
Panel B: 3x3 smoke-threshold x heat-percentile grid, spline adjustment,
         ZIP code-day clustering (read from results_import/rev2_matrix3x3.csv)."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)

# ---- Panel A data: (label, beta, {cov: (lo, hi)}, wcb)
A_ROWS = [
    ("No continuous adjustment", -1.78,
     {"person": (-3.13, -0.43), "zip3": (-3.13, -0.43), "zip3day": (-3.03, -0.53),
      "twoway": (-3.21, -0.35)}, "0.010"),
    ("Linear mean temperature", -1.49,
     {"person": (-2.84, -0.14), "zip3": (-2.86, -0.13), "zip3day": (-2.74, -0.24),
      "twoway": (-2.92, -0.06)}, "0.028"),
    ("Linear + quadratic", -1.48,
     {"person": (-2.83, -0.12), "zip3": (-2.84, -0.11), "zip3day": (-2.73, -0.23),
      "twoway": (-2.91, -0.05)}, "--"),
    ("Spline, mean temp. (primary)", -1.44,
     {"person": (-2.80, -0.09), "zip3": (-2.78, -0.10), "zip3day": (-2.69, -0.20),
      "twoway": (-2.87, -0.02)}, "0.028"),
    ("Splines, mean + overnight min.", -1.45,
     {"person": (-2.81, -0.10), "zip3": (-2.81, -0.09), "zip3day": (-2.70, -0.21),
      "twoway": (-2.88, -0.03)}, "--"),
    ("Spline, overnight min.", -1.54,
     {"person": (-2.89, -0.18), "zip3": (-2.93, -0.15), "zip3day": (-2.79, -0.28),
      "twoway": (-2.97, -0.10)}, "0.024"),
    ("Hot night on overnight min.", -0.79,
     {"person": (-1.91, 0.34), "zip3": (-2.05, 0.47), "zip3day": (-1.83, 0.26),
      "twoway": (-1.97, 0.40)}, "0.212"),
]
COVS = [("person", "Participant", "#4d4d4d", 0.27),
        ("zip3", "ZIP code", "#74add1", 0.09),
        ("zip3day", "ZIP code-day", "#b2182b", -0.09),
        ("twoway", "Two-way", "#8073ac", -0.27)]

plt.rcParams.update({"font.family": "Helvetica", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "axes.linewidth": 0.8})

fig = plt.figure(figsize=(12.8, 5.0))

# ================================================================ Panel A
axA = fig.add_axes([0.185, 0.13, 0.30, 0.72])
NA = len(A_ROWS)
ysA = list(range(NA))[::-1]
for y, (lab, b, cis, wcb) in zip(ysA, A_ROWS):
    for cov, _, color, off in COVS:
        lo, hi = cis[cov]
        lw = 1.7 if cov in ("person", "zip3day") else 1.1
        axA.plot([lo, hi], [y + off, y + off], color=color, lw=lw,
                 solid_capstyle="round", zorder=2)
        if cov == "person":
            axA.plot(b, y + off, "s", color=color, ms=5, zorder=3)
        elif cov == "zip3day":
            axA.plot(b, y + off, "o", color=color, ms=4, zorder=3)
axA.axvline(0, color="0.45", lw=0.8, ls="--", zorder=1)
axA.set_yticks(ysA)
axA.set_yticklabels([r[0] for r in A_ROWS], fontsize=8.6)
axA.tick_params(axis="y", length=0)
axA.set_xlabel("Smoke × heat interaction (minutes of sleep)", fontsize=8.5)
axA.set_xlim(-3.5, 0.8)
axA.set_ylim(-0.65, NA - 0.35)
for y in range(NA - 1):
    axA.axhline(y + 0.5, color="0.9", lw=0.6, zorder=0)

posA = axA.get_position()
def rowyA(y):
    return posA.y0 + (y - (-0.65)) / ((NA - 0.35) - (-0.65)) * posA.height
X_BETA, X_WCB = 0.497, 0.617
for y, (lab, b, cis, wcb) in zip(ysA, A_ROWS):
    lo, hi = cis["zip3day"]
    fig.text(X_BETA, rowyA(y), f"{b:.2f} ({lo:.2f}, {hi:.2f})", fontsize=8.2, va="center")
    sig = wcb not in ("--",) and float(wcb) < 0.05
    fig.text(X_WCB, rowyA(y), wcb, fontsize=8.2, va="center", ha="left",
             color="#b2182b" if sig else "0.25",
             fontweight="bold" if sig else "normal")
hdrA = posA.y1 + 0.045
fig.text(X_BETA, hdrA, "β (95% CI), ZIP code-day", fontsize=8.2, fontweight="bold")
fig.text(X_WCB, hdrA, "Wild bootstrap P", fontsize=8.2, fontweight="bold")
fig.text(0.02, 0.945, "A  Temperature specification and variance clustering",
         fontsize=11, fontweight="bold")

handles = [plt.Line2D([0], [0], color=c, lw=2) for _, _, c, _ in COVS]
fig.legend(handles, [n for _, n, _, _ in COVS], loc="lower left",
           bbox_to_anchor=(0.497, 0.015), frameon=False, fontsize=7.4,
           ncol=4, handlelength=1.3, columnspacing=1.0,
           title="Clustering level", title_fontsize=7.6)

# ================================================================ Panel B
m = pd.read_csv(f"{BASE}/results_import/rev2_matrix3x3.csv")
m["lo"] = m.int_beta - 1.96 * m.se_zip3day
m["hi"] = m.int_beta + 1.96 * m.se_zip3day
order = [(10, "P85"), (10, "P90"), (10, "P95"),
         (20, "P85"), (20, "P90"), (20, "P95"),
         (35, "P85"), (35, "P90"), (35, "P95")]
axB = fig.add_axes([0.755, 0.13, 0.225, 0.72])
NB = len(order)
ysB = []
y = NB + 1.6
for i, (thr, hq) in enumerate(order):
    if i % 3 == 0:
        y -= 0.9
    y -= 1.0
    ysB.append(y)
for (thr, hq), y in zip(order, ysB):
    r = m[(m.smoke_thr == thr) & (m.heat_q == hq)].iloc[0]
    sig = r.p_zip3day < 0.05
    color = "#b2182b" if sig else "#5a5a5a"
    axB.plot([r.lo, r.hi], [y, y], color=color, lw=1.7, solid_capstyle="round", zorder=2)
    axB.plot(r.int_beta, y, "o", color=color, ms=4.5, zorder=3)
    star = ("P90" if (thr, hq) == (20, "P90") else "")
    axB.text(r.hi + 0.25, y, f"{r.int_beta:.2f}", fontsize=7.6, va="center", color="0.25")
axB.axvline(0, color="0.45", lw=0.8, ls="--", zorder=1)
axB.set_yticks(ysB)
axB.set_yticklabels([f"{hq[1:]}th pct." for _, hq in order], fontsize=8.2)
axB.tick_params(axis="y", length=0)
axB.set_xlabel("Smoke × heat interaction (minutes of sleep)", fontsize=8.5)
axB.set_xlim(-9.6, 2.6)
axB.set_ylim(min(ysB) - 0.8, max(ysB) + 0.8)

posB = axB.get_position()
def rowyB(yv):
    lo, hi = min(ysB) - 0.8, max(ysB) + 0.8
    return posB.y0 + (yv - lo) / (hi - lo) * posB.height
for gi, thr in enumerate((10, 20, 35)):
    ytop = ysB[gi * 3]
    fig.text(posB.x0 - 0.085, rowyB(ytop) + 0.035,
             f"Smoke >{thr} µg/m³", fontsize=8.4, fontweight="bold")
fig.text(0.60, 0.945, "B  Threshold grid, spline adjustment (ZIP code-day CI)",
         fontsize=11, fontweight="bold")
fig.text(posB.x0 - 0.085, 0.015,
         "Red: P<0.05, ZIP code-day clustering. Primary cell: >20 µg/m³ × 90th percentile.",
         fontsize=6.8, color="0.35")

for ext, dpi in (("png", 600), ("svg", None), ("pdf", None)):
    fig.savefig(f"{BASE}/figures/Fig4_interaction_robustness.{ext}",
                dpi=dpi if dpi else "figure")
print("saved Fig4_interaction_robustness.*")
