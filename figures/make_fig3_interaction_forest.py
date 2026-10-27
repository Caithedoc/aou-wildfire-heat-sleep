#!/usr/bin/env python3
"""New Figure 3: prespecified subgroup interaction forest.
Point estimates/CIs/n from Phase E2 (published Table 4 values); formal
P_interaction from results_import/rev_subgroup_pint.csv."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)
pint = pd.read_csv(f"{BASE}/results_import/rev_subgroup_pint.csv")
pint = pint[pint.term.str.startswith("sx_")].set_index("model")

def pv(model):
    p = pint.loc[model, "p"]
    return "<0.001" if p < 0.001 else f"{p:.3f}"

ROWS = [
    ("Cardiorespiratory disease*", -1.48, -2.43, -0.54, "6,094",  pv("pint_cardioresp")),
    ("Obstructive sleep apnea",   -1.13, -2.18, -0.09, "5,437",  pv("pint_osa_f")),
    ("No obstructive sleep apnea",-0.33, -0.76,  0.11, "24,610", "—"),
    ("Age <65 years",             -0.65, -1.17, -0.14, "20,090", pv("pint_age65")),
    ("Age ≥65 years",        -0.11, -0.76,  0.54, "10,550", "—"),
    ("Most deprived ADI quartile",-0.32, -1.17,  0.54, "8,345",  pv("pint_deprived_top25")),
]
FULL = (-0.58, -1.02, -0.14)

plt.rcParams.update({"font.family": "Helvetica", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "axes.linewidth": 0.8})

fig = plt.figure(figsize=(9.2, 4.3))
ax = fig.add_axes([0.235, 0.155, 0.345, 0.685])

NY = len(ROWS)
ys = list(range(NY))[::-1]
for y, (lab, b, lo, hi, n, pi) in zip(ys, ROWS):
    sig = hi < 0
    color = "#b2182b" if sig else "#5a5a5a"
    ax.plot([lo, hi], [y, y], color=color, lw=1.6, solid_capstyle="round", zorder=2)
    ax.plot(b, y, "s", color=color, ms=6, zorder=3)

ax.axvline(0, color="0.45", lw=0.8, ls="--", zorder=1)
ax.axvspan(FULL[1], FULL[2], color="#f2b26b", alpha=0.18, zorder=0)
ax.axvline(FULL[0], color="#e08214", lw=1.1, ls=":", zorder=1)

ax.set_yticks(ys)
ax.set_yticklabels([r[0] for r in ROWS], fontsize=9)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("Change in total sleep time on smoke nights (minutes)", fontsize=8.5)
ax.set_xlim(-2.9, 1.3)
ax.set_ylim(-0.6, NY - 0.4)

pos = ax.get_position()
def rowy(y):
    return pos.y0 + (y - (-0.6)) / ((NY - 0.4) - (-0.6)) * pos.height

X_BETA, X_N, X_P = 0.615, 0.845, 0.885
for y, (lab, b, lo, hi, n, pi) in zip(ys, ROWS):
    fig.text(X_BETA, rowy(y), f"{b:.2f} ({lo:.2f}, {hi:.2f})", fontsize=8.5, va="center")
    fig.text(X_N, rowy(y), n, fontsize=8.5, va="center", ha="right")
    sigp = pi not in ("—",) and (pi == "<0.001" or float(pi) < 0.05)
    fig.text(X_P, rowy(y), pi, fontsize=8.5, va="center",
             color="#b2182b" if sigp else "0.25",
             fontweight="bold" if sigp else "normal")

hdr_y = pos.y1 + 0.05
fig.text(X_BETA, hdr_y, "β (95% CI), min", fontsize=8.5, fontweight="bold")
fig.text(X_N, hdr_y, "Participants", fontsize=8.5, fontweight="bold", ha="right")
fig.text(X_P, hdr_y, "P interaction", fontsize=8.5, fontweight="bold")

fig.text(0.045, 0.945, "Smoke-night association by subgroup",
         fontsize=11, fontweight="bold")
fig.text(0.235, 0.022,
         "Orange band and dotted line: full-sample estimate (-0.58 min; 95% CI: -1.02, -0.14).\n"
         "P interaction from smoke × subgroup product terms in the full sample; —, complement category of the contrast above.\nThe cardiorespiratory composite was defined post hoc; other subgroups were prespecified.",
         fontsize= 7, color="0.35")

for ext, dpi in (("png", 600), ("svg", None), ("pdf", None)):
    fig.savefig(f"{BASE}/figures/Fig3_interaction_forest.{ext}",
                dpi=dpi if dpi else "figure")
print("saved Fig3_interaction_forest.*")
