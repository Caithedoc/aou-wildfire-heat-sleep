"""17_min_nights_sweep: 每人最少分析夜数阈值扫描 (仿 Nat Health DST 文 Supp Fig 3), 2026-09-28.

在 AoU Jupyter VM 里跑: `nohup python -u 17_min_nights_sweep.py > rev6_sweep.log 2>&1 &`  (约 5-10 分钟)
读同目录缓存 analysis_nights.parquet; 设计变量与 12_recluster_tempflex.py 完全一致 (P90 按 zip3 x 月在分析夜上取,
smoke > 20, 样条结点 5/27.5/50/72.5/95 在全样本上定, 各阈值子样本沿用全样本结点)。
对每个阈值 k (每人分析夜数 >= k), 重拟合 M3 (sd + hd + sxh + tmean 样条) 与四格模型, 报 sxh 与 excess,
SE 用 zip3-日聚类 (正文报告口径) 与 person 聚类各一份。
输出:
  rev6_sweep_curve.csv     k = 1..400: 达到 >= k 夜的人数、夜数、复合夜数 (A 面板曲线)
  rev6_sweep_estimates.csv 阈值 k in K: 交互 sxh / excess 的 beta, se, lo, hi, p, 样本量 (B 面板)
1-20 计数置 NaN。"""
import time, numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy.stats import norm
T0 = time.time(); THRESH = 20; MIN_CELL = 20
K_CURVE = np.arange(1, 401); K_EST = [1, 5, 10, 20, 30, 60, 90, 120, 180, 270, 365]
def log(*a): print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)

COLS = ["person_id", "sleep_date", "ym", "dow", "zip3", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]: df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["date"].dt.month.astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
df["zip3_date"] = df["zip3"] + "_" + df["date"].dt.strftime("%Y%m%d")
df["sd"] = (df["smoke"] > THRESH).astype(float)
df["hd"] = (df["tmean"] >= df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
df["sxh"] = df["sd"] * df["hd"]
df["g_h"] = ((df.hd == 1) & (df.sd == 0)).astype(float); df["g_s"] = ((df.hd == 0) & (df.sd == 1)).astype(float); df["g_b"] = df["sxh"]
log("rows", len(df), "persons", df["person_id"].nunique())

def ncs_basis(x, k):
    x = np.asarray(x, float); k = np.asarray(k, float)
    def d(j): return (np.maximum(x - k[j], 0) ** 3 - np.maximum(x - k[-1], 0) ** 3) / (k[-1] - k[j])
    return np.column_stack([x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])
KN = np.percentile(df["tmean"], [5, 27.5, 50, 72.5, 95]); S = ncs_basis(df["tmean"].values, KN)
TS = [f"ts{i}" for i in range(S.shape[1])]
for i, nm in enumerate(TS): df[nm] = S[:, i]
log("knots", np.round(KN, 2))

def cluster_cov(X, u, codes, XtX_inv):
    G = codes.max() + 1; n, k = X.shape
    Sg = np.column_stack([np.bincount(codes, weights=X[:, j] * u, minlength=G) for j in range(k)])
    return (G / (G - 1)) * ((n - 1) / (n - k)) * (XtX_inv @ (Sg.T @ Sg) @ XtX_inv)

def fit(d, xs, cvec_map):
    d = d[d.groupby("stratum")["minute_asleep"].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d["minute_asleep"] - g["minute_asleep"].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    XtX_inv = np.linalg.inv(X.T @ X); beta = XtX_inv @ (X.T @ Y); u = Y - X @ beta
    covs = {"person": cluster_cov(X, u, pd.factorize(d["person_id"])[0], XtX_inv),
            "zip3_day": cluster_cov(X, u, pd.factorize(d["zip3_date"])[0], XtX_inv)}
    rows = []
    for nm, cvec in cvec_map.items():
        b = float(cvec @ beta)
        for cv, V in covs.items():
            se = float(np.sqrt(cvec @ V @ cvec))
            rows.append({"term": nm, "cov": cv, "beta": b, "se": se, "lo": b - 1.96 * se, "hi": b + 1.96 * se,
                         "p": 2 * norm.sf(abs(b / se))})
    return pd.DataFrame(rows), len(d), d["person_id"].nunique(), int(d["sxh"].sum())

# A: curve
pn = df.groupby("person_id").agg(n=("minute_asleep", "size"), both=("sxh", "sum"))
curve = pd.DataFrame({"k": K_CURVE,
                      "n_persons": [(pn["n"] >= k).sum() for k in K_CURVE],
                      "n_nights": [pn.loc[pn["n"] >= k, "n"].sum() for k in K_CURVE],
                      "n_both": [pn.loc[pn["n"] >= k, "both"].sum() for k in K_CURVE]})
for c in ["n_persons", "n_nights", "n_both"]: curve.loc[(curve[c] > 0) & (curve[c] < MIN_CELL), c] = np.nan
curve.to_csv("rev6_sweep_curve.csv", index=False); log("curve done")

# B: estimates
xs_int = ["sd", "hd", "sxh"] + TS; xs4 = ["g_h", "g_s", "g_b"] + TS
c_int = np.zeros(len(xs_int)); c_int[2] = 1
c_exc = np.zeros(len(xs4)); c_exc[:3] = (-1, -1, 1)
out = []
for k in K_EST:
    keep = pn.index[pn["n"] >= k]; d = df[df["person_id"].isin(keep)]
    if d["sxh"].sum() < MIN_CELL: log(f"k={k}: too few joint nights, skip"); continue
    r1, n_n, n_p, n_b = fit(d, xs_int, {"sxh": c_int}); r2, *_ = fit(d, xs4, {"excess": c_exc})
    r = pd.concat([r1, r2]); r["k"] = k; r["n_nights"] = n_n; r["n_persons"] = n_p; r["n_both"] = n_b
    out.append(r); log(f"k={k}: persons {n_p:,} nights {n_n:,} both {n_b:,} | sxh {r1.beta.iloc[0]:+.2f} (zip3day se {r1.se.iloc[1]:.2f})")
res = pd.concat(out)
for c in ["n_nights", "n_persons", "n_both"]: res.loc[(res[c] > 0) & (res[c] < MIN_CELL), c] = np.nan
res.to_csv("rev6_sweep_estimates.csv", index=False); log("done")
