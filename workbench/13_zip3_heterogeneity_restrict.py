"""13_zip3_heterogeneity_restrict: ZCTA->ZIP3 聚合误差的限定敏感性分析 (2026-09-22).

在 AoU Jupyter VM 里跑: `nohup python -u 13_zip3_heterogeneity_restrict.py > rev3.log 2>&1 &`
需要同目录: analysis_nights.parquet (缓存, 同 11/12) + zip3_heterogeneity.csv
(本地 local_exposure_prep/07_zip3_aggregation_error.py 产出, 纯公开数据派生, 上传即可)。

事前规定 (见 07 脚本头): 估计量 = 主模型 M3 (sd + hd + sxh + tmean 自然样条),
报 person / zip3 / zip3-日 / 双向 四种聚类 + ZIP3 wild cluster bootstrap。
热日阈值 (zip3 x 月 P90) 与样条节点均在全样本上定义, 限定只删 ZIP3, 不改暴露定义。
  FULL  全样本 (复现 rev2 M3, 作为核对: sxh 应为 -1.44)
  R1    H_smoke <= 人口加权中位数 (07 给出)
  R1n   H_smoke <= 分析夜加权中位数 (本脚本计算)
  R2    H_smoke 最低三分位
  R3    剔除 H_smoke 最高三分位
  EM    全样本上 sxh x H_smoke(标准化) 效应修饰
输出 (全部聚合, n<20 单元置 NaN):
  rev3_restrict.csv      各限定的交互项 / 四格 excess, 四种聚类
  rev3_wcb.csv           各限定 sxh 的 ZIP3 WCB p
  rev3_em.csv            效应修饰模型
  rev3_exposure_by_H.csv 按 H 三分位的夜数 / 复合夜数
  rev3_meta.csv
"""

import time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

THRESH = 20
MIN_CELL = 20
RNG = 20260922
B_WCB = 999
T0 = time.time()


def log(*a):
    print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


# ---------------------------------------------------------------- 读缓存 (同 12)
COLS = ["person_id", "sleep_date", "ym", "dow", "zip3", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["month"] = df["date"].dt.month
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                 + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["month"].astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
df["zip3_date"] = df["zip3"] + "_" + df["date"].dt.strftime("%Y%m%d")
log("rows:", len(df), "| persons:", df["person_id"].nunique(), "| zip3:", df["zip3"].nunique())

het = pd.read_csv("zip3_heterogeneity.csv", dtype={"zip3": str})
het["zip3"] = het["zip3"].str.zfill(3)
df = df.merge(het[["zip3", "H_smoke", "H_compound", "R1_le_popwt_median",
                   "R2_bottom_tercile", "R3_drop_top_tercile"]], on="zip3", how="left")
cov_H = df["H_smoke"].notna().mean()
log(f"H_smoke 覆盖率 (夜) = {cov_H:.4f}")
df = df[df["H_smoke"].notna()].copy()


# ---------------------------------------------------------------- 估计机件 (同 12)
def ncs_basis(x, knots):
    x = np.asarray(x, float); k = np.asarray(knots, float)
    def d(j):
        return ((np.maximum(x - k[j], 0) ** 3 - np.maximum(x - k[-1], 0) ** 3)
                / (k[-1] - k[j]))
    return np.column_stack([x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])


def cluster_cov(X, u, codes, XtX_inv):
    G = codes.max() + 1
    n, k = X.shape
    S = np.column_stack([np.bincount(codes, weights=X[:, j] * u, minlength=G)
                         for j in range(k)])
    return (G / (G - 1)) * ((n - 1) / (n - k)) * (XtX_inv @ (S.T @ S) @ XtX_inv)


def hc1_cov(X, u, XtX_inv):
    n, k = X.shape
    return (n / (n - k)) * (XtX_inv @ ((X * (u ** 2)[:, None]).T @ X) @ XtX_inv)


def fit_multi(d, y, xs, label, contrasts=None, wcb_term=None):
    d = d.dropna(subset=[y] + xs).copy()
    d = d[d.groupby("stratum")[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ Y)
    u = Y - X @ beta
    cp = pd.factorize(d["person_id"])[0]
    cz = pd.factorize(d["zip3"])[0]
    czd = pd.factorize(d["zip3_date"])[0]
    Vp, Vz, Vzd = (cluster_cov(X, u, c, XtX_inv) for c in (cp, cz, czd))
    covs = {"person": Vp, "zip3": Vz, "zip3_day": Vzd,
            "twoway_p_z3d": Vp + Vzd - hc1_cov(X, u, XtX_inv)}
    items = [(nm, np.eye(len(xs))[i]) for i, nm in enumerate(xs)]
    if contrasts:
        items += list(contrasts.items())
    rows = []
    for cvname, V in covs.items():
        for nm, cvec in items:
            b = float(cvec @ beta)
            var = float(cvec @ V @ cvec)
            se = np.sqrt(var) if var > 0 else np.nan
            rows.append({"model": label, "term": nm, "cov": cvname, "beta": b, "se": se,
                         "lo": b - 1.96 * se, "hi": b + 1.96 * se,
                         "p": 2 * norm.sf(abs(b / se)) if se == se and se > 0 else np.nan,
                         "n_nights": len(d), "n_persons": d["person_id"].nunique(),
                         "n_clusters_zip3": int(cz.max() + 1),
                         "n_clusters_zip3day": int(czd.max() + 1)})
    wcb = None
    if wcb_term is not None:
        j = xs.index(wcb_term)
        se_z = np.sqrt(Vz[j, j]); t_obs = beta[j] / se_z
        Xr = X[:, [i for i in range(len(xs)) if i != j]]
        br = np.linalg.lstsq(Xr, Y, rcond=None)[0]
        ur, fit0 = Y - Xr @ br, Xr @ br
        rng = np.random.default_rng(RNG)
        G = cz.max() + 1
        AtY = XtX_inv @ X.T
        exceed = 0
        for _ in range(B_WCB):
            yb = fit0 + ur * rng.choice([-1.0, 1.0], size=G)[cz]
            bb = AtY @ yb
            Vb = cluster_cov(X, yb - X @ bb, cz, XtX_inv)
            exceed += abs(bb[j] / np.sqrt(Vb[j, j])) >= abs(t_obs)
        wcb = {"model": label, "term": wcb_term, "B": B_WCB, "t_obs_zip3": t_obs,
               "p_wcb": (exceed + 1) / (B_WCB + 1), "beta": beta[j], "se_zip3": se_z}
        log(f"  WCB {label}: t={t_obs:.3f}, p_wcb={wcb['p_wcb']:.4f}")
    return pd.DataFrame(rows), wcb


# ---------------------------------------------------------------- 暴露 (全样本定义)
df["sd"] = (df["smoke"] > THRESH).astype(float)
df["hd"] = (df["tmean"] >= df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
df["sxh"] = df["sd"] * df["hd"]
df["g_h"] = ((df.hd == 1) & (df.sd == 0)).astype(float)
df["g_s"] = ((df.hd == 0) & (df.sd == 1)).astype(float)
df["g_b"] = ((df.hd == 1) & (df.sd == 1)).astype(float)
kn = np.percentile(df["tmean"], [5, 27.5, 50, 72.5, 95])
Sb = ncs_basis(df["tmean"].values, kn)
TS = [f"ts{i}" for i in range(Sb.shape[1])]
for i, nm in enumerate(TS):
    df[nm] = Sb[:, i]

# 分析夜加权中位数 (限有烟雾日的 ZIP3, 与 07 的人口加权口径平行)
ever = df.groupby("zip3")["sd"].transform("max") > 0
h = df.loc[ever, "H_smoke"].sort_values().values
med_n = float(h[len(h) // 2])
df["R1n"] = (df["H_smoke"] <= med_n).astype(int)
log(f"H_smoke 夜加权中位数 = {med_n:.4f}")

# ---------------------------------------------------------------- 按 H 三分位的暴露分布
df["H_tercile"] = np.select([df["R2_bottom_tercile"] == 1, df["R3_drop_top_tercile"] == 1],
                            ["T1_low", "T2_mid"], "T3_high")
ex = df.groupby("H_tercile").agg(nights=("sd", "size"), smoke_nights=("sd", "sum"),
                                 compound_nights=("g_b", "sum"), zip3=("zip3", "nunique"),
                                 persons=("person_id", "nunique"))
for c in ex.columns:
    ex.loc[ex[c] < MIN_CELL, c] = np.nan
ex.to_csv("rev3_exposure_by_H.csv")
log("=== 按 H 三分位 ===\n", ex.to_string())

# ---------------------------------------------------------------- 限定分析
SETS = [("FULL", np.ones(len(df), bool)),
        ("R1_popwt_median", df["R1_le_popwt_median"].values == 1),
        ("R1n_nightwt_median", df["R1n"].values == 1),
        ("R2_bottom_tercile", df["R2_bottom_tercile"].values == 1),
        ("R3_drop_top_tercile", df["R3_drop_top_tercile"].values == 1)]
res, wcbs = [], []
for lab, m in SETS:
    d = df[m]
    n_b = int(d["g_b"].sum())
    if n_b < MIN_CELL:
        log(f"{lab}: 复合夜 {n_b} < {MIN_CELL}, 跳过"); continue
    oi, w = fit_multi(d, "minute_asleep", ["sd", "hd", "sxh"] + TS, f"int_{lab}", wcb_term="sxh")
    cvec = np.zeros(3 + len(TS)); cvec[:3] = [-1, -1, 1]
    o4, _ = fit_multi(d, "minute_asleep", ["g_h", "g_s", "g_b"] + TS, f"four_{lab}",
                      contrasts={"excess_b_minus_h_minus_s": cvec})
    for o in (oi, o4):
        o["subset"], o["n_both"] = lab, n_b
        res.append(o)
    w["subset"], w["n_both"] = lab, n_b
    wcbs.append(w)
    k = oi[oi.term == "sxh"].set_index("cov")
    log(f"{lab}: n_both={n_b} sxh={k.loc['person','beta']:+.3f} "
        f"p person={k.loc['person','p']:.4f} zip3day={k.loc['zip3_day','p']:.4f} "
        f"zip3={k.loc['zip3','p']:.4f} twoway={k.loc['twoway_p_z3d','p']:.4f}")
pd.concat(res).to_csv("rev3_restrict.csv", index=False)
pd.DataFrame(wcbs).to_csv("rev3_wcb.csv", index=False)

# ---------------------------------------------------------------- 效应修饰 (连续 H)
df["Hz"] = (df["H_smoke"] - df["H_smoke"].mean()) / df["H_smoke"].std()
for c in ("sd", "hd", "sxh"):
    df[f"{c}_Hz"] = df[c] * df["Hz"]
em, _ = fit_multi(df, "minute_asleep", ["sd", "hd", "sxh", "sd_Hz", "hd_Hz", "sxh_Hz"] + TS, "em_Hz")
em.to_csv("rev3_em.csv", index=False)
log("=== EM (person) ===\n", em[em["cov"] == "person"][["term", "beta", "lo", "hi", "p"]]
    .head(6).round(4).to_string(index=False))

pd.DataFrame([{"key": "rows", "value": len(df)}, {"key": "H_night_coverage", "value": cov_H},
              {"key": "H_nightwt_median", "value": med_n},
              {"key": "tmean_knots", "value": ";".join(f"{x:.2f}" for x in kn)},
              {"key": "runtime_s", "value": int(time.time() - T0)}]).to_csv("rev3_meta.csv", index=False)
log("ALL DONE — 导出前人工过一遍 rev3_*.csv 确认无 1-20 计数单元。")
