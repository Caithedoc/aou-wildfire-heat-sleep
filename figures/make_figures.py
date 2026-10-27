"""野火烟 x 睡眠 主图集 (2026-08-28)。
风格: Okabe-Ito 色盲安全, Helvetica/Arial 无衬线, forest 三栏 (标签|图区|数值)。
输出: PNG 600dpi + SVG。数据 = Workbench 导出的聚合系数 (无个体级数据)。
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).parent
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})
# Okabe-Ito
TEAL, SKY, VERM, BLUE, ORNG, GREY = "#009E73", "#56B4E9", "#D55E00", "#0072B2", "#E69F00", "#666666"


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def forest(ax, labels, beta, lo, hi, colors=None, ref=0.0):
    y = np.arange(len(labels))[::-1]
    colors = colors or [BLUE] * len(labels)
    ax.axvline(ref, color=GREY, lw=0.8, ls="--", zorder=1)
    for i, (b, l, h, c) in enumerate(zip(beta, lo, hi, colors)):
        ax.plot([l, h], [y[i], y[i]], color=c, lw=1.6, solid_capstyle="round", zorder=2)
        ax.plot([b], [y[i]], "o", color=c, ms=5.5, zorder=3,
                mec="white", mew=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0)
    return y


# ---------------------------------------------------------------- Fig 1 稳健性
rows = [
    ("Primary (smoke >20, lag 0-1)",              -0.580, -1.023, -0.136, VERM),
    ("Without temperature adjustment",            -0.756, -1.199, -0.312, GREY),
    ("Wear-probability weighted (IPW)",           -0.438, -0.883,  0.007, BLUE),
    ("Threshold >10 ug/m3",                       -0.068, -0.427,  0.291, GREY),
    ("Threshold >35 ug/m3",                       -0.491, -1.252,  0.271, GREY),
    ("Wider sleep window (1-16 h)",               -0.203, -0.656,  0.250, ORNG),
    ("Excluding hot nights",                      -0.350, -0.783,  0.083, ORNG),
    ("Excluding 2023",                            -0.546, -1.177,  0.084, GREY),
    ("Western US (ZIP3 8xx-9xx)",                 -0.050, -0.820,  0.720, GREY),
    ("Eastern US (ZIP3 0xx-2xx)",                 -0.577, -1.430,  0.280, GREY),
    ("June 2023 East Coast episode",              -0.144, -1.230,  0.940, GREY),
    ("Negative control: next-day smoke",           0.144, -0.306,  0.595, TEAL),
    ("Negative control: low-smoke years",          1.063, -0.423,  2.549, TEAL),
]
rows = [r for r in rows if not np.isnan(r[1])]
labels = [r[0] for r in rows]
fig, ax = plt.subplots(figsize=(6.6, 4.2))
forest(ax, labels, [r[1] for r in rows], [r[2] for r in rows],
       [r[3] for r in rows], [r[4] for r in rows])
ax.set_xlabel("Change in total sleep time on smoke nights (minutes)")
ax.set_xlim(-1.6, 2.7)
for i, r in enumerate(rows):
    ax.text(2.75, len(rows) - 1 - i, f"{r[1]:.2f} ({r[2]:.2f}, {r[3]:.2f})",
            va="center", ha="left", fontsize=7.5, color="#222222")
ax.set_title("Smoke-night association is small and not robust across specifications",
             loc="left", fontsize=9.5, pad=8)
save(fig, "Fig1_robustness_forest")

# ---------------------------------------------------------------- Fig 2 复合暴露 (头条)
fig, ax = plt.subplots(figsize=(5.2, 3.4))
cells = ["Neither", "Hot only", "Smoky only", "Hot + smoky"]
vals = [0.0, -1.718, -0.569, -4.021]
los = [0.0, -2.121, -0.992, -5.260]
his = [0.0, -1.316, -0.145, -2.782]
cols = [GREY, ORNG, SKY, VERM]
x = np.arange(4)
ax.axhline(0, color=GREY, lw=0.8)
ax.bar(x, vals, color=cols, width=0.6, zorder=2)
for i in range(1, 4):
    ax.plot([x[i], x[i]], [los[i], his[i]], color="#222222", lw=1.2, zorder=3)
add = vals[1] + vals[2]
ax.plot([x[3] - 0.42, x[3] + 0.42], [add, add], color="#222222", lw=1.4, ls=":", zorder=4)
ax.annotate("additive\nexpectation", xy=(x[3] + 0.44, add), fontsize=7.2,
            va="center", ha="left", color="#222222")
ax.annotate("", xy=(x[3] - 0.30, add), xytext=(x[3] - 0.30, vals[3]),
            arrowprops=dict(arrowstyle="<->", lw=1.0, color="#222222"))
ax.text(x[3] - 0.36, (add + vals[3]) / 2, "excess\n-1.7 min", fontsize=7.2,
        ha="right", va="center", color="#222222")
ax.set_xticks(x)
ax.set_xticklabels(cells)
ax.set_ylabel("Change in total sleep time (minutes)")
ax.set_ylim(-5.8, 0.7)
ax.set_title("Compound heat-smoke nights: effect exceeds the sum of its parts",
             loc="left", fontsize=9.5, pad=8)
ax.text(0.02, -5.55, "Interaction term -1.44 (95% CI -2.80, -0.09), p = 0.037   |   "
        "18,168 compound nights", fontsize=7, color="#444444")
save(fig, "Fig2_compound_exposure")

# ---------------------------------------------------------------- Fig 3 滞后曲线
fig, ax = plt.subplots(figsize=(5.4, 3.0))
lag = np.arange(8)
b = [-0.432, 0.222, 0.022, -0.285, 0.584, 0.732, -0.008, 0.331]
lo = [-0.890, -0.262, -0.460, -0.765, 0.097, 0.244, -0.502, -0.135]
hi = [0.026, 0.705, 0.503, 0.194, 1.071, 1.221, 0.486, 0.797]
ax.axhline(0, color=GREY, lw=0.8, ls="--")
ax.fill_between(lag, lo, hi, color=SKY, alpha=0.25, lw=0)
ax.plot(lag, b, "-o", color=BLUE, lw=1.5, ms=4.5, mec="white", mew=0.8)
ax.set_xlabel("Days since smoke exposure (lag)")
ax.set_ylabel("Change in total sleep time (minutes)")
ax.set_title("No coherent recovery trajectory after smoke exposure",
             loc="left", fontsize=9.5, pad=8)
save(fig, "Fig3_distributed_lag")

# ---------------------------------------------------------------- Fig 4 亚组
subs = [
    ("Cardiorespiratory disease",  -1.483, -2.425, -0.541, 6094,  VERM),
    ("Obstructive sleep apnoea",   -1.134, -2.181, -0.088, 5437,  VERM),
    ("Age <65 y",                  -0.653, -1.165, -0.140, 20090, BLUE),
    ("Age >=65 y",                 -0.109, -0.757,  0.538, 10550, GREY),
    ("No sleep apnoea",            -0.326, -0.762,  0.110, 24610, GREY),
    ("Most deprived quartile",     -0.318, -1.172,  0.536, 8345,  GREY),
]
fig, ax = plt.subplots(figsize=(6.2, 2.7))
forest(ax, [s[0] for s in subs], [s[1] for s in subs], [s[2] for s in subs],
       [s[3] for s in subs], [s[5] for s in subs])
ax.set_xlabel("Change in total sleep time on smoke nights (minutes)")
ax.set_xlim(-2.7, 1.5)
for i, s in enumerate(subs):
    ax.text(1.52, len(subs) - 1 - i, f"{s[1]:.2f} ({s[2]:.2f}, {s[3]:.2f})   n={s[4]:,}",
            va="center", ha="left", fontsize=7.2, color="#222222")
ax.set_title("Cardiorespiratory disease marks the clearest susceptibility",
             loc="left", fontsize=9.5, pad=8)
save(fig, "Fig4_subgroups")

# ---------------------------------------------------------------- Fig 5 GATES 失败
fig, ax = plt.subplots(figsize=(4.8, 3.2))
q = np.arange(1, 6)
pred = [-8.189, -1.904, -0.841, 0.402, 6.250]
obs = [-2.643, -0.652, -1.089, 0.260, 1.587]
olo = [-3.613, -1.458, -1.856, -0.537, 0.630]
ohi = [-1.672, 0.153, -0.321, 1.058, 2.543]
ax.axhline(0, color=GREY, lw=0.8, ls="--")
ax.plot(q, pred, "-s", color=GREY, lw=1.2, ms=4.5, label="Predicted by causal forest")
ax.plot(q, obs, "-o", color=VERM, lw=1.6, ms=5, mec="white", mew=0.8,
        label="Observed (GATES)")
for i in range(5):
    ax.plot([q[i], q[i]], [olo[i], ohi[i]], color=VERM, lw=1.2)
ax.set_xticks(q)
ax.set_xlabel("Quintile of predicted individual effect")
ax.set_ylabel("Change in total sleep time (minutes)")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title("Causal-forest heterogeneity fails calibration:\n"
             "the top quintile implies smoke *increases* sleep",
             loc="left", fontsize=9, pad=8)
save(fig, "Fig5_gates_calibration")

print("\nall figures written to", OUT)
