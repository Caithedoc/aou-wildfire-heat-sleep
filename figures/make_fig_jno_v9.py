#!/usr/bin/env python3
"""JNO V9 figures (2026-09-22). JAMA style: no in-plot titles, italic capital P without leading zero,
true minus signs, "3-digit ZIP code" wording. Outputs PNG (600 dpi) + SVG + PDF.

Figure 1 (new): A, 4-category joint-exposure estimates vs neither (Table 3, person-clustered CIs);
                B, smoke x heat interaction across the 3x3 threshold grid (eTable 10 panel B,
                spline temperature, ZIP code-day CIs, matching the primary interaction row).
Figure 2 (was the single "Figure"): continuous joint smoke x temperature response curves
                (results_import/rev_joint_curves.csv), restyled; same data as V8.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from lu_palette import seq, ACC, SEQ, BAND_ALPHA

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (release edit; original was the author's local project path)
OUT = f"{BASE}/figures"
INK, INK2, GRID = "#1a1a1a", "#555555", "#d9d9d9"
RAMP = seq(4)  # Lu 2025 palette: viridis samples for ordered temperature percentiles

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 8.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.7, "axes.edgecolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
    "axes.unicode_minus": True, "svg.fonttype": "none", "pdf.fonttype": 42,
})


def m(x, nd=2):
    """Number with a true minus sign."""
    return f"{x:.{nd}f}".replace("-", "−")


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=600, bbox_inches="tight",
                    metadata={"Creator": "Zhang, Cai"} if ext != "png" else None)
    plt.close(fig)
    print("saved", name)


# ---------------------------------------------------------------- Figure 1
CATS = [  # label, nights, beta, lo, hi   (Table 3)
    ("Hot only", 176_696, -1.69, -2.09, -1.29),
    ("Smoky only", 174_282, -0.57, -0.99, -0.14),
    ("Hot and smoky", 17_295, -4.03, -5.27, -2.79),
]
ADDITIVE, EXCESS = -2.25, (-1.78, -3.13, -0.43)
GRID_B = {  # (smoke thr, heat pct): beta, lo, hi   (eTable 10 panel B, ZIP code-day CI)
    (10, 85): (-0.78, -1.49, -0.07), (10, 90): (-0.68, -1.52, 0.15), (10, 95): (-0.54, -1.65, 0.56),
    (20, 85): (-1.51, -2.56, -0.46), (20, 90): (-1.44, -2.69, -0.20), (20, 95): (-2.06, -3.62, -0.51),
    (35, 85): (-4.72, -6.49, -2.95), (35, 90): (-3.92, -5.88, -1.96), (35, 95): (-6.15, -8.67, -3.63),
}

fig = plt.figure(figsize=(7.2, 3.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1], wspace=0.42)

# --- A: point ranges
ax = fig.add_subplot(gs[0])
ys = np.arange(len(CATS))[::-1]
for y, (lab, n, b, lo, hi) in zip(ys, CATS):
    joint = lab == "Hot and smoky"
    col = ACC if joint else INK2
    ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=3)
    ax.scatter([b], [y], s=34, color=col, zorder=4, edgecolor="white", linewidth=1.0)
    ax.text(0.35, y, f"{m(b)} ({m(lo)} to {m(hi)})", va="center", ha="left", fontsize=7.5, color=INK)
yj = ys[2]
ax.scatter([ADDITIVE], [yj], s=34, marker="D", facecolor="white", edgecolor=INK2, linewidth=1.0, zorder=4)
ax.annotate("", xy=(-4.03, yj - 0.32), xytext=(ADDITIVE, yj - 0.32),
            arrowprops=dict(arrowstyle="<->", color=ACC, lw=0.9, shrinkA=0, shrinkB=0))
ax.text((ADDITIVE - 4.03) / 2, yj - 0.52,
        f"Excess beyond additivity {m(EXCESS[0])}\n(95% CI, {m(EXCESS[1])} to {m(EXCESS[2])})",
        ha="center", va="top", fontsize=7, color=ACC)
ax.text(ADDITIVE, yj + 0.24, "Additive\nexpectation", ha="center", va="bottom", fontsize=6.8, color=INK2)
ax.axvline(0, color=INK2, lw=0.7, ls="--", zorder=1)
ax.set_yticks(ys)
ax.set_yticklabels([f"{lab}\n{n:,} nights" for lab, n, *_ in CATS], fontsize=8)
ax.set_xlim(-6.0, 3.9)
ax.set_xticks([-6, -4, -2, 0])
ax.set_ylim(-1.3, len(CATS) - 0.4)
ax.set_xlabel("Difference in total sleep time vs nights\nwith neither exposure, min")
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_bounds(-6, 0)
ax.tick_params(axis="y", length=0)
ax.text(0.35, len(CATS) - 0.55, "Estimate (95% CI)", fontsize=7.5, color=INK, fontweight="bold")
ax.text(-0.02, 1.02, "A", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")

# --- B: 3x3 heatmap
bx = fig.add_subplot(gs[1])
smk, heat = [10, 20, 35], [85, 90, 95]
Z = np.array([[GRID_B[(s, h)][0] for h in heat] for s in smk])
im = bx.imshow(-Z, cmap=SEQ, vmin=0, vmax=6.5, origin="lower", aspect="auto")
for i, s in enumerate(smk):
    for j, h in enumerate(heat):
        b, lo, hi = GRID_B[(s, h)]
        dark = -b < 4.2  # viridis: dark purple/teal below, yellow above
        bx.text(j, i + 0.1, m(b), ha="center", va="center", fontsize=9, fontweight="bold",
                color="white" if dark else INK)
        bx.text(j, i - 0.2, f"{m(lo)} to {m(hi)}", ha="center", va="center", fontsize=6.3,
                color="white" if dark else INK2)
bx.add_patch(Rectangle((0.5, 0.5), 1, 1, fill=False, edgecolor="white", lw=1.6, zorder=5))
bx.set_xticks(range(3)); bx.set_xticklabels([f"{h}th" for h in heat])
bx.set_yticks(range(3)); bx.set_yticklabels([f">{s}" for s in smk])
bx.set_xlabel("Hot-night threshold, percentile")
bx.set_ylabel("Smoke-night threshold, µg/m³")
for sp in bx.spines.values():
    sp.set_visible(False)
bx.tick_params(length=0)
cb = fig.colorbar(im, ax=bx, fraction=0.05, pad=0.03)
cb.set_ticks([0, 2, 4, 6]); cb.set_ticklabels(["0", "−2", "−4", "−6"])
cb.set_label("Smoke × heat interaction, min", fontsize=7.5)
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, length=0)
bx.text(-0.02, 1.02, "B", transform=bx.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")
save(fig, "Fig_JNO1_joint_categories_grid")

# ---------------------------------------------------------------- Figure 2
df = pd.read_csv(f"{BASE}/results_import/rev_joint_curves.csv")
SPEC = [(0.50, "50th (median night)"), (0.75, "75th"), (0.90, "90th (hot-night threshold)"), (0.95, "95th")]
fig, ax = plt.subplots(figsize=(7.0, 4.3))
for (tp, lab), col in zip(SPEC, RAMP):
    d = df[df.tpct == tp]
    ax.fill_between(d.smoke, d.lo, d.hi, color=col, alpha=0.12, lw=0, zorder=1)
for (tp, lab), col in zip(SPEC, RAMP):
    d = df[df.tpct == tp]
    ax.plot(d.smoke, d.dTST, color=col, lw=2.0, label=lab, zorder=3)
    last = d.iloc[-1]
    ax.text(last.smoke + 0.8, last.dTST, f"{int(tp * 100)}th", color=INK, fontsize=7.5, va="center")
ax.axhline(0, color=INK2, lw=0.7, ls="--", zorder=0)
ax.axvline(20, color="#bdbdbd", lw=0.7, ls=":", zorder=0)
ax.text(20.6, -3.25, "Primary smoke-night\nthreshold (20 µg/m³)", fontsize=7, color=INK2, va="bottom")
ax.set_xlim(0, 64)
ax.set_xticks(range(0, 61, 10))
ax.set_xlabel("Wildfire-specific PM$_{2.5}$, µg/m³")
ax.set_ylabel("Change in total sleep time vs no smoke, min")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
leg = ax.legend(title="Temperature percentile within\n3-digit ZIP code and calendar month",
                frameon=False, fontsize=7.5, title_fontsize=7.5, loc="lower left")
leg._legend_box.align = "left"
ax.text(0.99, 0.97, "Smoke × temperature interaction, $\\it{P}$ = .01\n(spline joint model)",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK2)
save(fig, "Fig_JNO2_joint_response_curves")
