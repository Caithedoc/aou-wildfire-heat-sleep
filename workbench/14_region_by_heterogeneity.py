"""14_region_by_heterogeneity: rev3 最低异质性三分位零效应的地区分解 (2026-09-22, 事后分析).

背景: rev3 (13_zip3_heterogeneity_restrict.py) 中只留 H_smoke 最低三分位时交互项 +0.15
(CI -2.29, 2.58), 而剔除最高三分位时 -1.72。本地 07 显示低 H 的 ZIP3 集中在东部。
问题: 最低三分位的零效应是异质性 (暴露误差) 造成的, 还是地区构成 (东部远程输送烟雾) 造成的?
**本分析在看到 rev3 结果后设计, 属事后分析, 文中须如此标注。**

在 AoU Jupyter VM 里跑: `nohup python -u 14_region_by_heterogeneity.py > rev4.log 2>&1 &`
需要: analysis_nights.parquet + zip3_heterogeneity.csv (同 13)。
地区 = ZIP 首位 0-2 East / 3-7 Central / 8-9 West (同 11_revision.py:308)。
输出 (聚合, n<20 置 NaN):
  rev4_composition.csv   地区 x H 三分位 的夜 / 烟雾夜 / 复合夜 / ZIP3 / 人数
  rev4_region.csv        A 各地区交互项 (四种聚类) + B 地区 x (T1 vs T2+T3)
  rev4_wcb.csv           A 各地区 WCB
  rev4_C_em_*.csv        C 全样本: sxh x H 修饰, 同时控制 sxh x 地区
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


f = df["zip3"].str[0].astype(int)
df["region"] = np.select([f <= 2, f <= 7], ["East", "Central"], "West")
df["H_tercile"] = np.select([df["R2_bottom_tercile"] == 1, df["R3_drop_top_tercile"] == 1],
                            ["T1_low", "T2_mid"], "T3_high")

# ---------------------------------------------------------------- 构成
comp = df.groupby(["region", "H_tercile"]).agg(
    nights=("sd", "size"), smoke_nights=("sd", "sum"), compound_nights=("g_b", "sum"),
    zip3=("zip3", "nunique"), persons=("person_id", "nunique"))
for c in comp.columns:
    comp.loc[comp[c] < MIN_CELL, c] = np.nan
comp.to_csv("rev4_composition.csv")
log("=== 地区 x H 三分位 ===\n", comp.to_string())

XS = ["sd", "hd", "sxh"] + TS
res, wcbs = [], []


def run(d, lab, wcb):
    n_b = int(d["g_b"].sum())
    if n_b < MIN_CELL:
        log(f"{lab}: 复合夜 < {MIN_CELL}, 跳过"); return
    o, w = fit_multi(d, "minute_asleep", XS, lab, wcb_term="sxh" if wcb else None)
    o["n_both"] = n_b
    res.append(o)
    if w:
        w["n_both"] = n_b; wcbs.append(w)
    k = o[o.term == "sxh"].set_index("cov")
    log(f"{lab}: n_both={n_b} sxh={k.loc['zip3','beta']:+.3f} "
        f"CI_zip3=({k.loc['zip3','lo']:+.2f},{k.loc['zip3','hi']:+.2f}) "
        f"p person={k.loc['person','p']:.3f} zip3={k.loc['zip3','p']:.3f} "
        f"zip3day={k.loc['zip3_day','p']:.3f}")


# A: 各地区
for rg in ("East", "Central", "West"):
    run(df[df.region == rg], f"A_{rg}", wcb=True)
# B: 地区 x (T1 vs T2+T3)
for rg in ("East", "Central", "West"):
    d = df[df.region == rg]
    run(d[d.H_tercile == "T1_low"], f"B_{rg}_T1", wcb=False)
    run(d[d.H_tercile != "T1_low"], f"B_{rg}_T23", wcb=False)
pd.concat(res).to_csv("rev4_region.csv", index=False)
pd.DataFrame(wcbs).to_csv("rev4_wcb.csv", index=False)

# C: H 修饰, 控制地区修饰 (East 为参照)
df["Hz"] = (df["H_smoke"] - df["H_smoke"].mean()) / df["H_smoke"].std()
extra = []
for c in ("sd", "hd", "sxh"):
    df[f"{c}_Hz"] = df[c] * df["Hz"]
    for rg in ("Central", "West"):
        df[f"{c}_{rg}"] = df[c] * (df["region"] == rg)
        extra.append(f"{c}_{rg}")
for lab, xs in (("C_em_H_only", ["sd", "hd", "sxh", "sd_Hz", "hd_Hz", "sxh_Hz"] + TS),
                ("C_em_region_only", ["sd", "hd", "sxh"] + extra + TS),
                ("C_em_H_plus_region", ["sd", "hd", "sxh", "sd_Hz", "hd_Hz", "sxh_Hz"] + extra + TS)):
    o, _ = fit_multi(df, "minute_asleep", xs, lab)
    res_c = o[o["term"].str.startswith("sxh")]
    extra_rows = res_c[res_c["cov"].isin(["person", "zip3"])]
    log(f"=== {lab} ===\n", extra_rows[["term", "cov", "beta", "lo", "hi", "p"]].round(4).to_string(index=False))
    o.to_csv(f"rev4_{lab}.csv", index=False)
log("ALL DONE — 导出前人工检查 rev4_*.csv 无 1-20 计数单元。")
