"""20_sleep_domains_joint: 四域睡眠结局 joint 模型 (2026-10-12, ANALYSIS_PLAN_ADDENDUM_SLEEPDOMAINS).

在 AoU Jupyter VM 里跑 (先跑 19 重建缓存):
    nohup python -u 20_sleep_domains_joint.py > rev8_domains.log 2>&1 &

模型机件与 18_secondary_joint.py 逐行相同 (正文 M3 规格); 只扩展结局与夜集:
  时长:   minute_asleep                       (时长集 = 全部分析夜; 锚点 excess 应复现 -1.78)
  连续性: eff_pp, minute_awake, n_awake5      (前两者时长集; n_awake5 结构集)
  时相:   onset_min, offset_min, midsleep_min, midsleep_shift, minute_to_fall_asleep (结构集)
  结构:   minute_light/deep/rem 分钟 + pct_light/deep/rem 占比 (结构集)
  复合:   healthy_night, sh_score, healthy_night_s (时长集), healthy_night_v2 (结构集)
结构集 = 分期分钟非空的夜; 每个结局的 N 单独报 (Liao 2025 两套 N 的做法)。
缺列的结局自动跳过并写进 run_meta (不会静默为零: 启动时打印 OUTS 实际清单)。

输出: rev8_domains_coefs.csv / rev8_domains_wcb.csv / rev8_domains_cells.csv / rev8_domains_bh.csv / rev8_run_meta.csv
"""
import os
import time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

T0 = time.time()
THRESH = 20
MIN_CELL = 20
RNG = 20261012
B_WCB = int(os.environ.get("B_WCB", 999))

DOMAIN = {  # 结局 -> (域, 夜集)
    "minute_asleep": ("duration", "all"),
    "eff_pp": ("continuity", "all"), "minute_awake": ("continuity", "all"), "n_awake5": ("continuity", "stage"),
    "onset_min": ("timing", "stage"), "offset_min": ("timing", "stage"), "midsleep_min": ("timing", "stage"),
    "midsleep_shift": ("timing", "stage"), "minute_to_fall_asleep": ("timing", "stage"),
    "minute_light": ("architecture", "stage"), "minute_deep": ("architecture", "stage"), "minute_rem": ("architecture", "stage"),
    "pct_light": ("architecture", "stage"), "pct_deep": ("architecture", "stage"), "pct_rem": ("architecture", "stage"),
    "healthy_night": ("composite", "all"), "sh_score": ("composite", "all"), "healthy_night_s": ("composite", "all"),
    "healthy_night_v2": ("composite", "stage"),
}


def log(*a):
    print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


# ---------------------------------------------------------------- 读缓存 (列按实际存在取)
t = pq.read_table("analysis_nights.parquet")
BASE = ["person_id", "sleep_date", "ym", "dow", "zip3", "smoke", "smoke_l1", "tmean"]
RAW = ["minute_asleep", "minute_awake", "minute_in_bed", "minute_to_fall_asleep", "minute_restless",
       "minute_light", "minute_deep", "minute_rem", "minute_wake", "sleep_efficiency",
       "onset_min", "offset_min", "midsleep_min", "midsleep_shift_raw", "gap_prev", "n_awake5", "n_awake_any", "n_awake_all"]
miss_base = [c for c in BASE if c not in t.column_names]
assert not miss_base, f"parquet 缺基础列: {miss_base}"
have_raw = [c for c in RAW if c in t.column_names]
log("parquet 可用结局列:", have_raw, "| 缺:", [c for c in RAW if c not in t.column_names])
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in BASE + have_raw})
for c in have_raw + ["smoke", "smoke_l1", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")

# ---------------------------------------------------------------- 导出变量 (附录 §1 定义)
df["eff_pp"] = df["sleep_efficiency"] * 100 if df["sleep_efficiency"].max() <= 1.5 else df["sleep_efficiency"]
df["mode_stages"] = df["minute_wake"].notna().astype(float) if "minute_wake" in df else np.nan
stage_ok = df["minute_deep"].notna() & df["minute_rem"].notna() & df["minute_light"].notna() if "minute_light" in df \
    else df["minute_deep"].notna() & df["minute_rem"].notna()
for s in ["light", "deep", "rem"]:
    if f"minute_{s}" in df:
        df[f"pct_{s}"] = 100.0 * df[f"minute_{s}"] / df["minute_asleep"]
        df.loc[~stage_ok | (df["minute_asleep"] <= 0), f"pct_{s}"] = np.nan
if "midsleep_shift_raw" in df:
    df["midsleep_shift"] = df["midsleep_shift_raw"].where(df["gap_prev"] == 1)
# 复合 v1 (SLEEPHEALTH 附录, 与 18 脚本逐行相同)
df["c_dur"] = ((df["minute_asleep"] >= 420) & (df["minute_asleep"] <= 540)).astype(float)
df["c_eff"] = (df["eff_pp"] >= 85).astype(float)
df["c_waso"] = (df["minute_awake"] <= 20).astype(float)
df["sh_score"] = df["c_dur"] + df["c_eff"] + df["c_waso"]
df["healthy_night"] = 100.0 * (df["sh_score"] == 3)
_miss = df[["minute_asleep", "eff_pp", "minute_awake"]].isna().any(axis=1)
df.loc[_miss, ["c_dur", "c_eff", "c_waso", "sh_score", "healthy_night"]] = np.nan
df["healthy_night_s"] = 100.0 * (df["c_dur"].eq(1) & (df["eff_pp"] >= 80) & (df["minute_awake"] <= 30))
df.loc[_miss, "healthy_night_s"] = np.nan
# 复合 v2: v1 + 潜伏期 <= 30 + n_awake5 <= 1 (结构集)
if "minute_to_fall_asleep" in df and "n_awake5" in df:
    v2_ok = ~_miss & df["minute_to_fall_asleep"].notna() & df["n_awake5"].notna()
    df["healthy_night_v2"] = 100.0 * ((df["sh_score"] == 3) & (df["minute_to_fall_asleep"] <= 30) & (df["n_awake5"] <= 1))
    df.loc[~v2_ok, "healthy_night_v2"] = np.nan

# 结构集旗标: 分期分钟非空
df["in_stage_set"] = stage_ok.astype(float)
for y, (dom, nset) in DOMAIN.items():
    if y in df and nset == "stage":
        df.loc[~stage_ok, y] = np.nan      # 结构集结局只在分期夜上定义

OUTS = [y for y in DOMAIN if y in df.columns and df[y].notna().sum() > 0]
skipped = [y for y in DOMAIN if y not in OUTS]
log("OUTS 实际运行:", OUTS); log("跳过 (缺列或全空):", skipped)

df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["date"].dt.month.astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
df["zip3_date"] = df["zip3"] + "_" + df["date"].dt.strftime("%Y%m%d")
log("rows", len(df), "| persons", df["person_id"].nunique(), "| stage-set nights", int(df["in_stage_set"].sum()))

# ---------------------------------------------------------------- 设计变量 (同 12/18)
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


KN = np.percentile(df["tmean"], [5, 27.5, 50, 72.5, 95])     # 全样本结点, 结构集沿用 (附录 §3)
S = ncs_basis(df["tmean"].values, KN)
TS = [f"ts{i}" for i in range(S.shape[1])]
for i, nm in enumerate(TS):
    df[nm] = S[:, i]
log("tmean knots", np.round(KN, 2))


# ---------------------------------------------------------------- 估计机件 (同 18)
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
    if len(d) < 1000 or d["sxh"].sum() < MIN_CELL:
        log(f"  skip {label}/{y}: nights={len(d)} both={int(d['sxh'].sum())}")
        return pd.DataFrame(), None
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
            rows.append({"model": label, "outcome": y, "domain": DOMAIN[y][0], "night_set": DOMAIN[y][1],
                         "term": nm, "cov": cv, "beta": b, "se": se,
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
        wcb = {"model": label, "outcome": y, "domain": DOMAIN[y][0], "term": wcb_term, "B": B_WCB,
               "beta": beta[j], "t_obs_zip3": t_obs, "p_wcb": (exceed + 1) / (B_WCB + 1)}
        log(f"  WCB {y}: p_wcb={wcb['p_wcb']:.4f}")
    return pd.DataFrame(rows), wcb


# ---------------------------------------------------------------- 拟合
coef, wcbs = [], []
x4 = ["g_h", "g_s", "g_b"] + TS
c_exc = np.zeros(len(x4)); c_exc[:3] = (-1, -1, 1)
for y in OUTS:
    n_y = int(df[y].notna().sum())
    log(f"=== {y} [{DOMAIN[y][0]} / {DOMAIN[y][1]}]: {n_y:,} nights ===")
    o, _ = fit_multi(df, y, x4, "four_M3", contrasts={"excess_b_minus_h_minus_s": c_exc}); coef.append(o)
    o, w = fit_multi(df, y, ["sd", "hd", "sxh"] + TS, "int_M3", wcb_term="sxh"); coef.append(o)
    if w: wcbs.append(w)
    o, _ = fit_multi(df, y, ["sd", "smoke_l1_day"] + TS, "marg_M3"); coef.append(o)
    o, _ = fit_multi(df, y, ["sd", "tmean"], "tlin"); coef.append(o)
    k = pd.concat(coef[-4:])
    if len(k):
        k = k[(k["cov"] == "zip3_day") & k["term"].isin(["g_h", "g_s", "g_b", "excess_b_minus_h_minus_s", "sxh", "tmean"])]
        log("\n" + k[["model", "term", "beta", "lo", "hi", "p", "n_nights"]].round(4).to_string(index=False))

coef = pd.concat([c for c in coef if len(c)], ignore_index=True)
for c in ["n_nights", "n_persons", "n_both"]:
    coef.loc[(coef[c] > 0) & (coef[c] < MIN_CELL), c] = np.nan
coef.to_csv("rev8_domains_coefs.csv", index=False)
pd.DataFrame(wcbs).to_csv("rev8_domains_wcb.csv", index=False)

# ---------------------------------------------------------------- 按格描述 (含夜集与设备模式)
cell = np.select([df.g_b == 1, df.g_h == 1, df.g_s == 1], ["both", "hot_only", "smoke_only"], default="neither")
rows = []
for y in OUTS + ["mode_stages", "in_stage_set"]:
    for cname, s in df.assign(cell=cell).groupby("cell")[y]:
        n = int(s.notna().sum())
        rows.append({"outcome": y, "cell": cname, "n_nights": n if n >= MIN_CELL else np.nan,
                     "mean": s.mean() if n >= MIN_CELL else np.nan,
                     "sd": s.std() if n >= MIN_CELL else np.nan})
pd.DataFrame(rows).to_csv("rev8_domains_cells.csv", index=False)

# ---------------------------------------------------------------- BH 参考, 按域分族 (附录 §4.3)
sel = coef[(coef["cov"] == "zip3_day") & (coef["outcome"] != "minute_asleep")
           & (((coef["model"] == "four_M3") & (coef["term"] == "excess_b_minus_h_minus_s"))
              | ((coef["model"] == "int_M3") & (coef["term"] == "sxh")))].copy()
for (fam, dom), g in sel.groupby(["model", "domain"]):
    p = g["p"].values; o = np.argsort(p); m = len(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    sel.loc[g.index[o], "q_bh_domain"] = np.minimum(q, 1)
sel.to_csv("rev8_domains_bh.csv", index=False)

meta = [{"key": "rows", "value": len(df)}, {"key": "persons", "value": df["person_id"].nunique()},
        {"key": "n_both", "value": int(df["g_b"].sum())}, {"key": "stage_set_nights", "value": int(df["in_stage_set"].sum())},
        {"key": "stage_set_persons", "value": int(df.loc[df["in_stage_set"] == 1, "person_id"].nunique())},
        {"key": "mode_stages_pct", "value": round(100 * float(df["mode_stages"].mean()), 2) if "minute_wake" in df else "na"},
        {"key": "outs_run", "value": ";".join(OUTS)}, {"key": "outs_skipped", "value": ";".join(skipped)},
        {"key": "wcb_B", "value": B_WCB}, {"key": "tmean_knots", "value": ";".join(f"{x:.2f}" for x in KN)},
        {"key": "healthy_night_pct", "value": round(float(df["healthy_night"].mean()), 2)},
        {"key": "runtime_s", "value": int(time.time() - T0)}]
if "onset_min" in df:
    meta.append({"key": "onset_min_median", "value": round(float(df["onset_min"].median()), 1)})
pd.DataFrame(meta).to_csv("rev8_run_meta.csv", index=False)
log("ALL DONE. 锚点: minute_asleep four_M3 excess (zip3_day) 应约 -1.78 (-3.13, -0.43)。"
    " 导出前人工确认 rev8_*.csv 无 1-20 计数。")
