"""18_secondary_joint: 次要睡眠结局升格进正文 (2026-10-10, 回应 Qian 批注 108/109/158/185/241).

在 AoU Jupyter VM 里跑: `nohup python -u 18_secondary_joint.py > rev7_secondary.log 2>&1 &`
读同目录缓存 analysis_nights.parquet。设计变量与 12_recluster_tempflex.py 完全一致:
smoke > 20, 热 = tmean >= zip3 x 月 P90 (在全分析样本上取, 分期子样本沿用), tmean 自然样条 5 结点
(全样本 5/27.5/50/72.5/95 分位, 分期子样本沿用全样本结点)。主结局 TST 一并重跑作锚点,
excess 应复现正文 -1.78。

结局 (全部预设为 secondary, ANALYSIS_PLAN.md:15/91): WASO (minute_awake), 深睡 (minute_deep),
REM (minute_rem), 睡眠效率 (sleep_efficiency x 100 = 百分点)。
2026-10-12 追加 (ANALYSIS_PLAN_ADDENDUM_SLEEPHEALTH_20261012.md, 开跑前预设): healthy_night (NSF 三成分
全达标, 0/100 线性概率, 系数 = 百分点) 与 sh_score (0-3 达标项数)。

每个结局拟合三套 (均含 tmean 样条, 即正文 M3 规格):
  four   : g_h + g_s + g_b              -> 热单独 / 烟单独 / 联合 + excess (b - h - s)
  int    : sd + hd + sxh                -> 交互项 sxh, 外加 zip3 层级 WCB p
  marg   : smoke_day + smoke_l1 + 样条  -> 边际烟 (复刻 Table 2 次级行, 换成样条温度)
  tlin   : sd + tmean 线性              -> 每 1 °C 温度主效应 (Qian 254: 热主效应表)
SE: person / zip3 / zip3-日 / 双向 四种; 正文报 zip3-日 (与主结局口径一致)。

输出 (全部聚合, 1-20 计数置 NaN):
  rev7_secondary_coefs.csv   model|outcome|term|cov|beta|se|lo|hi|p|n_nights|n_persons
  rev7_secondary_wcb.csv     各结局 sxh 的 WCB p
  rev7_secondary_cells.csv   各结局按四格的均值/SD/夜数 (Table 1 补充 + 效应量换算)
  rev7_secondary_bh.csv      zip3-日 口径下 excess 与 sxh 的 BH q 值 (仅作参考, 正文不做正式校正)
  rev7_run_meta.csv
"""
import time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

T0 = time.time()
THRESH = 20
MIN_CELL = 20
RNG = 20261010
B_WCB = int(__import__("os").environ.get("B_WCB", 999))
OUTS = ["minute_asleep", "minute_awake", "minute_deep", "minute_rem", "eff_pp",
        "healthy_night", "sh_score", "healthy_night_s"]   # 后三个: ADDENDUM_SLEEPHEALTH 2026-10-12


def log(*a):
    print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


# ---------------------------------------------------------------- 读缓存
COLS = ["person_id", "sleep_date", "ym", "dow", "zip3",
        "minute_asleep", "minute_awake", "minute_deep", "minute_rem", "sleep_efficiency",
        "smoke", "smoke_l1", "tmean"]
t = pq.read_table("analysis_nights.parquet")
miss = [c for c in COLS if c not in t.column_names]
assert not miss, f"parquet 缺列: {miss}"
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in COLS[5:]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["eff_pp"] = df["sleep_efficiency"] * 100 if df["sleep_efficiency"].max() <= 1.5 else df["sleep_efficiency"]

# ---------------------------------------------------------------- 预设复合结局 (ANALYSIS_PLAN_ADDENDUM_SLEEPHEALTH_20261012.md)
# NSF 共识切点: 时长 7-9 h (Hirshkowitz 2015) / 效率 >= 85% / WASO <= 20 min (Ohayon 2017)。
# v1 三成分; 潜伏期与觉醒次数待 04b 重建后在 20 脚本作 v2。任一成分缺失 -> 复合结局缺失。
df["c_dur"] = ((df["minute_asleep"] >= 420) & (df["minute_asleep"] <= 540)).astype(float)
df["c_eff"] = (df["eff_pp"] >= 85).astype(float)
df["c_waso"] = (df["minute_awake"] <= 20).astype(float)
df["sh_score"] = df["c_dur"] + df["c_eff"] + df["c_waso"]          # 0-3 达标项数
df["healthy_night"] = 100.0 * (df["sh_score"] == 3)                 # 0/100, 系数 = 百分点
_miss = df[["minute_asleep", "eff_pp", "minute_awake"]].isna().any(axis=1)
df.loc[_miss, ["c_dur", "c_eff", "c_waso", "sh_score", "healthy_night"]] = np.nan
# 预设敏感性: 设备偏倚放宽一档 (WASO <= 30, 效率 >= 80), 附录用, 不替换 v1
df["healthy_night_s"] = 100.0 * (df["c_dur"].eq(1) & (df["eff_pp"] >= 80) & (df["minute_awake"] <= 30))
df.loc[_miss, "healthy_night_s"] = np.nan
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                 + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["date"].dt.month.astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
df["zip3_date"] = df["zip3"] + "_" + df["date"].dt.strftime("%Y%m%d")
log("rows", len(df), "| persons", df["person_id"].nunique())

# ---------------------------------------------------------------- 设计变量 (同 12)
df["sd"] = (df["smoke"] > THRESH).astype(float)
df["smoke_l1_day"] = (df["smoke_l1"] > THRESH).astype(float)
df["hd"] = (df["tmean"] >= df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
df["sxh"] = df["sd"] * df["hd"]
df["g_h"] = ((df.hd == 1) & (df.sd == 0)).astype(float)
df["g_s"] = ((df.hd == 0) & (df.sd == 1)).astype(float)
df["g_b"] = df["sxh"]


def ncs_basis(x, k):
    x = np.asarray(x, float); k = np.asarray(k, float)
    def d(j):
        return (np.maximum(x - k[j], 0) ** 3 - np.maximum(x - k[-1], 0) ** 3) / (k[-1] - k[j])
    return np.column_stack([x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])


KN = np.percentile(df["tmean"], [5, 27.5, 50, 72.5, 95])
S = ncs_basis(df["tmean"].values, KN)
TS = [f"ts{i}" for i in range(S.shape[1])]
for i, nm in enumerate(TS):
    df[nm] = S[:, i]
log("tmean knots", np.round(KN, 2))


# ---------------------------------------------------------------- 估计机件 (同 12)
def cluster_cov(X, u, codes, XtX_inv):
    G = codes.max() + 1; n, k = X.shape
    Sg = np.column_stack([np.bincount(codes, weights=X[:, j] * u, minlength=G) for j in range(k)])
    return (G / (G - 1)) * ((n - 1) / (n - k)) * (XtX_inv @ (Sg.T @ Sg) @ XtX_inv)


def hc1_cov(X, u, XtX_inv):
    n, k = X.shape
    return (n / (n - k)) * (XtX_inv @ ((X * (u ** 2)[:, None]).T @ X) @ XtX_inv)


def fit_multi(d, y, xs, label, contrasts=None, wcb_term=None):
    d = d.dropna(subset=[y] + xs)
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
    for cv, V in covs.items():
        for nm, cvec in items:
            b = float(cvec @ beta); var = float(cvec @ V @ cvec)
            se = np.sqrt(var) if var > 0 else np.nan
            rows.append({"model": label, "outcome": y, "term": nm, "cov": cv, "beta": b, "se": se,
                         "lo": b - 1.96 * se, "hi": b + 1.96 * se,
                         "p": 2 * norm.sf(abs(b / se)) if se == se and se > 0 else np.nan,
                         "n_nights": len(d), "n_persons": d["person_id"].nunique(),
                         "n_both": int(d["sxh"].sum())})
    wcb = None
    if wcb_term is not None:
        j = xs.index(wcb_term)
        t_obs = beta[j] / np.sqrt(Vz[j, j])
        keep = [i for i in range(len(xs)) if i != j]
        Xr = X[:, keep]
        br = np.linalg.lstsq(Xr, Y, rcond=None)[0]
        fit0 = Xr @ br; ur = Y - fit0
        AtY = XtX_inv @ X.T
        rng = np.random.default_rng(RNG); G = cz.max() + 1; exceed = 0
        for _ in range(B_WCB):
            yb = fit0 + ur * rng.choice([-1.0, 1.0], size=G)[cz]
            bb = AtY @ yb
            Vb = cluster_cov(X, yb - X @ bb, cz, XtX_inv)
            exceed += abs(bb[j] / np.sqrt(Vb[j, j])) >= abs(t_obs)
        wcb = {"model": label, "outcome": y, "term": wcb_term, "B": B_WCB,
               "beta": beta[j], "t_obs_zip3": t_obs, "p_wcb": (exceed + 1) / (B_WCB + 1)}
        log(f"  WCB {y}: p_wcb={wcb['p_wcb']:.4f}")
    return pd.DataFrame(rows), wcb


# ---------------------------------------------------------------- 拟合
coef, wcbs = [], []
x4 = ["g_h", "g_s", "g_b"] + TS
c_exc = np.zeros(len(x4)); c_exc[:3] = (-1, -1, 1)
for y in OUTS:
    n_y = int(df[y].notna().sum())
    log(f"=== {y}: {n_y:,} nights ===")
    o, _ = fit_multi(df, y, x4, "four_M3", contrasts={"excess_b_minus_h_minus_s": c_exc})
    coef.append(o)
    o, w = fit_multi(df, y, ["sd", "hd", "sxh"] + TS, "int_M3", wcb_term="sxh")
    coef.append(o); wcbs.append(w)
    o, _ = fit_multi(df, y, ["sd", "smoke_l1_day"] + TS, "marg_M3")
    coef.append(o)
    o, _ = fit_multi(df, y, ["sd", "tmean"], "tlin")
    coef.append(o)
    k = pd.concat(coef[-4:])
    k = k[(k["cov"] == "zip3_day") & k["term"].isin(["g_h", "g_s", "g_b", "excess_b_minus_h_minus_s", "sxh", "tmean"])]
    log("\n" + k[["model", "term", "beta", "lo", "hi", "p"]].round(4).to_string(index=False))

coef = pd.concat(coef, ignore_index=True)
for c in ["n_nights", "n_persons", "n_both"]:
    coef.loc[(coef[c] > 0) & (coef[c] < MIN_CELL), c] = np.nan
coef.to_csv("rev7_secondary_coefs.csv", index=False)
pd.DataFrame(wcbs).to_csv("rev7_secondary_wcb.csv", index=False)

# ---------------------------------------------------------------- 按格描述
cell = np.select([df.g_b == 1, df.g_h == 1, df.g_s == 1], ["both", "hot_only", "smoke_only"], default="neither")
rows = []
for y in OUTS:
    for cname, s in df.assign(cell=cell).groupby("cell")[y]:
        n = int(s.notna().sum())
        rows.append({"outcome": y, "cell": cname, "n_nights": n if n >= MIN_CELL else np.nan,
                     "mean": s.mean() if n >= MIN_CELL else np.nan,
                     "sd": s.std() if n >= MIN_CELL else np.nan})
pd.DataFrame(rows).to_csv("rev7_secondary_cells.csv", index=False)

# ---------------------------------------------------------------- BH 参考 (仅次要结局)
sel = coef[(coef["cov"] == "zip3_day") & (coef["outcome"] != "minute_asleep")
           & (((coef["model"] == "four_M3") & (coef["term"] == "excess_b_minus_h_minus_s"))
              | ((coef["model"] == "int_M3") & (coef["term"] == "sxh")))].copy()
for fam, g in sel.groupby("model"):
    p = g["p"].values; o = np.argsort(p); m = len(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    sel.loc[g.index[o], "q_bh"] = np.minimum(q, 1)
sel.to_csv("rev7_secondary_bh.csv", index=False)

pd.DataFrame([{"key": "rows", "value": len(df)},
              {"key": "persons", "value": df["person_id"].nunique()},
              {"key": "n_both", "value": int(df["g_b"].sum())},
              {"key": "wcb_B", "value": B_WCB},
              {"key": "tmean_knots", "value": ";".join(f"{x:.2f}" for x in KN)},
              {"key": "healthy_night_pct", "value": round(float(df["healthy_night"].mean()), 2)},
              {"key": "sh_score_mean", "value": round(float(df["sh_score"].mean()), 3)},
              {"key": "healthy_night_s_pct", "value": round(float(df["healthy_night_s"].mean()), 2)},
              {"key": "c_dur_pct", "value": round(100 * float(df["c_dur"].mean()), 2)},
              {"key": "c_eff_pct", "value": round(100 * float(df["c_eff"].mean()), 2)},
              {"key": "c_waso_pct", "value": round(100 * float(df["c_waso"].mean()), 2)},
              {"key": "runtime_s", "value": int(time.time() - T0)}]).to_csv("rev7_run_meta.csv", index=False)
log("ALL DONE. 锚点检查: minute_asleep four_M3 excess (zip3_day) 应约 -1.78 (-3.13, -0.43)。"
    " 导出前人工确认 rev7_*.csv 无 1-20 计数。")
