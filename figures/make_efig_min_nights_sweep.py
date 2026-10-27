"""eFigure: sensitivity of the smoke x heat interaction to the minimum number of analytic nights per participant,
styled after Nat Health DST paper Supplementary Fig 3 (A: sample-size curve with threshold lines; B: estimate bands).
Run: /usr/bin/python3 figures/make_efig_min_nights_sweep.py [curve.csv estimates.csv]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 9})
cf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results_import", "rev6_sweep_curve.csv")
ef = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results_import", "rev6_sweep_estimates.csv")
curve = pd.read_csv(cf); est = pd.read_csv(ef)
K = sorted(est["k"].unique()); colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(K)))
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), gridspec_kw={"width_ratios": [1.15, 1, 1]})
ax = axes[0]
ax.plot(curve["k"], curve["n_persons"] / 1000, color="#1f4e79", lw=1.4, label="Participants (thousands)")
ax2 = ax.twinx(); ax2.plot(curve["k"], curve["n_both"] / 1000, color="#c0504d", lw=1.4, ls="--", label="Hot-and-smoky nights (thousands)")
for k, c in zip(K, colors): ax.axvline(k, color=c, lw=0.8, ls=":")
ax.set_xlabel("Minimum analytic nights per participant, k"); ax.set_ylabel("Participants with ≥ k nights, thousands", color="#1f4e79")
ax2.set_ylabel("Hot-and-smoky nights retained, thousands", color="#c0504d"); ax.set_title("A. Sample retained", loc="left", fontsize=9)
ax.set_xlim(1, curve["k"].max())
for j, (term, ttl) in enumerate([("sxh", "B. Smoke × heat interaction (spline model)"), ("excess", "C. Excess beyond additivity (4-category model)")]):
    ax = axes[j + 1]; sub = est[(est["term"] == term) & (est["cov"] == "zip3_day")].sort_values("k")
    x = np.arange(len(sub))
    ax.fill_between(x, sub["lo"], sub["hi"], color="#bfbfbf", alpha=0.6, step="mid", label="95% CI (ZIP code-day clusters)")
    ax.plot(x, sub["beta"], color="black", lw=1.2, marker="o", ms=3.5)
    ax.axhline(0, color="black", lw=0.7, ls="--")
    ax.set_xticks(x); ax.set_xticklabels([f"{int(k)}\n{int(n):,}" if pd.notna(n) else str(int(k)) for k, n in zip(sub["k"], sub["n_persons"])], fontsize=7)
    ax.set_xlabel("Minimum analytic nights per participant, k (participants retained)"); ax.set_ylabel("Difference in total sleep time, min")
    ax.set_title(ttl, loc="left", fontsize=9)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "eFig_min_nights_sweep" + ("_PROXY" if "scratch" in cf else ""))
for ext in ("png", "svg", "pdf"): fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print("saved", out)
