"""12_recluster_tempflex: 编辑 panel 决定性再分析 (2026-08-30, EDITOR_PANEL findings A/B/C).

在 AoU Jupyter VM 里跑: `nohup python 12_recluster_tempflex.py > rev2.log 2>&1 &`
读 analysis_nights.parquet 缓存 (同 11_revision.py); 可选读 tmin_zip3_daily.csv.gz
(本地 05_tmin_zip3_prep.py 产出, 先上传到 VM 工作目录; 没有该文件时自动跳过 tmin 模型)。

回答三个问题:
  A. 聚类层级: 暴露在 zip3-日 赋值, 原 SE 只按 person 聚类 (Moulton)。
     每个模型同时报 person / zip3 / zip3-日 / 双向(person+zip3日, CGM) 四种 SE,
     外加 zip3 层级 wild cluster bootstrap (Rademacher, 限制式) p 值。
  B. 温度调整: 头条四分类模型原本无温度项。跑温度阶梯
     M0 无温度(复刻) / M1 tmean 线性 / M2 +二次 / M3 tmean 自然样条(5节点)
     / M4 tmean样条+tmin样条 / M5 tmin样条; 并导出按暴露格的温度分布表。
     tmin_night 用 PRISM 次日文件 (PRISM 日=UTC12止, D+1 文件覆盖 D->D+1 夜间最低温)。
  C. abstract 断言补模型: compound x cardioresp (及 OSA) 交互 + 分层四分类。

附带: 3x3 稳健性矩阵重跑(样条温度+新聚类)、compound 格的次级结局、
连续联合模型正主效应复查、超加性 excess 的 delta 法 CI。

输出 (全部聚合系数/计数, n<20 单元格置 NaN):
  rev2_run_meta.csv            运行元数据/样本量
  rev2_cell_temperature.csv    B: 按格温度分布 (决定 finding vs artifact 的那张表)
  rev2_fourcell_ladder.csv     B+A: 四分类系数+excess, 模型阶梯 x 聚类方式
  rev2_interaction_ladder.csv  B+A: 乘积项参数化, 模型阶梯 x 聚类方式
  rev2_wcb.csv                 A: wild cluster bootstrap p 值
  rev2_matrix3x3.csv           3x3 网格重跑
  rev2_subgroup_compound.csv   C: compound x cardioresp / osa
  rev2_secondary_compound.csv  compound 格次级结局
  rev2_joint_linear_check.csv  正烟雾主效应复查
"""

import os, subprocess, time
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
THRESH = 20
MIN_CELL = 20
RNG = 20260830
B_WCB = 999
T0 = time.time()

def log(*a):
    print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)

def bqcsv(sql, max_rows="100000000"):
    try:
        r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                            f"--max_rows={max_rows}", sql], capture_output=True, text=True)
    except Exception as e:
        print("BQ UNAVAILABLE:", e); return None
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[:500]); return None
    return pd.read_csv(StringIO(r.stdout))

# ---------------------------------------------------------------- 读缓存 (同 11)
COLS = ["person_id", "sleep_date", "ym", "dow", "src_id", "zip3",
        "minute_asleep", "minute_awake", "sleep_efficiency",
        "smoke", "smoke_l1", "tmean"]
t = pq.read_table("analysis_nights.parquet")
have = [c for c in COLS if c in t.column_names]
print("parquet cols missing:", [c for c in COLS if c not in t.column_names])
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in have})
for c in ["minute_asleep", "minute_awake", "sleep_efficiency", "smoke", "smoke_l1", "tmean"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["month"] = df["date"].dt.month
df["year"] = df["date"].dt.year
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["smoke_day"] = (df["smoke"] > THRESH).astype(float)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                 + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["month"].astype(str)
sz = df.groupby("stratum")["minute_asleep"].transform("size")
df = df[sz >= 2].copy()
df["zip3_date"] = df["zip3"] + "_" + df["date"].dt.strftime("%Y%m%d")
log("regression-consistent rows:", len(df), "| persons:", df["person_id"].nunique(),
    "| zip3:", df["zip3"].nunique(), "| zip3-days:", df["zip3_date"].nunique())

# ---------------------------------------------------------------- tmin (可选)
HAVE_TMIN = False
for f in ("tmin_zip3_daily.csv.gz", "tmin_zip3_daily.csv", "tmin_zip3_daily.parquet"):
    if os.path.exists(f):
        tm = (pd.read_parquet(f) if f.endswith("parquet") else
              pd.read_csv(f, dtype={"zip3": str}))
        tm["zip3"] = tm["zip3"].astype(str).str.zfill(3)
        tm["date"] = pd.to_datetime(tm["date"])
        # PRISM 日界=UTC12 (~当地凌晨4-8点): D+1 文件的 tmin 覆盖 D 夜间; 同日文件作敏感性
        nxt = tm.copy(); nxt["date"] = nxt["date"] - pd.Timedelta(days=1)
        nxt = nxt.rename(columns={"tmin_popwt": "tmin_night"})
        df = df.merge(nxt[["zip3", "date", "tmin_night"]], on=["zip3", "date"], how="left")
        df = df.merge(tm.rename(columns={"tmin_popwt": "tmin_same"})[
            ["zip3", "date", "tmin_same"]], on=["zip3", "date"], how="left")
        HAVE_TMIN = df["tmin_night"].notna().mean() > 0.95
        log(f"tmin merged from {f}: night coverage "
            f"{df['tmin_night'].notna().mean():.4f}, same-day {df['tmin_same'].notna().mean():.4f}"
            f" -> HAVE_TMIN={HAVE_TMIN}")
        break
if not HAVE_TMIN:
    log("tmin 不可用 — M4/M5 及 tmin 相关输出将跳过")

# ---------------------------------------------------------------- 估计机件
def ncs_basis(x, knots):
    x = np.asarray(x, float); k = np.asarray(knots, float)
    def d(j):
        return ((np.maximum(x - k[j], 0) ** 3 - np.maximum(x - k[-1], 0) ** 3)
                / (k[-1] - k[j]))
    return np.column_stack([x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])

def add_spline(d, col, prefix):
    kn = np.percentile(d[col].dropna(), [5, 27.5, 50, 72.5, 95])
    S = ncs_basis(d[col].values, kn)
    names = [f"{prefix}{i}" for i in range(S.shape[1])]
    for i, nm in enumerate(names):
        d[nm] = S[:, i]
    return names, kn

def cluster_cov(X, u, codes, XtX_inv):
    G = codes.max() + 1
    n, k = X.shape
    S = np.column_stack([np.bincount(codes, weights=X[:, j] * u, minlength=G)
                         for j in range(k)])
    corr = (G / (G - 1)) * ((n - 1) / (n - k))
    return corr * (XtX_inv @ (S.T @ S) @ XtX_inv)

def hc1_cov(X, u, XtX_inv):
    n, k = X.shape
    meat = (X * (u ** 2)[:, None]).T @ X
    return (n / (n - k)) * (XtX_inv @ meat @ XtX_inv)

def fit_multi(d, y, xs, label, contrasts=None, wcb_term=None):
    """stratum 内去均值 OLS; SE 按 person / zip3 / zip3-日 / 双向 分别报。
    contrasts: {name: 系数向量} 线性组合 (delta 法)。wcb_term: 对该列做 WCB。
    返回 (系数行 DataFrame, wcb 行 dict 或 None)。"""
    d = d.dropna(subset=[y] + xs).copy()
    g = d.groupby("stratum")
    d = d[g[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX)
    beta = XtX_inv @ (X.T @ Y)
    u = Y - X @ beta
    cp = pd.factorize(d["person_id"])[0]
    cz = pd.factorize(d["zip3"])[0]
    czd = pd.factorize(d["zip3_date"])[0]
    Vp = cluster_cov(X, u, cp, XtX_inv)
    Vz = cluster_cov(X, u, cz, XtX_inv)
    Vzd = cluster_cov(X, u, czd, XtX_inv)
    V2 = Vp + Vzd - hc1_cov(X, u, XtX_inv)          # CGM 双向
    covs = {"person": Vp, "zip3": Vz, "zip3_day": Vzd, "twoway_p_z3d": V2}
    from scipy.stats import norm
    rows = []
    items = [(nm, np.eye(len(xs))[i]) for i, nm in enumerate(xs)]
    if contrasts:
        items += list(contrasts.items())
    for cvname, V in covs.items():
        for nm, cvec in items:
            b = float(cvec @ beta)
            var = float(cvec @ V @ cvec)
            se = np.sqrt(var) if var > 0 else np.nan
            rows.append({"model": label, "outcome": y, "term": nm, "cov": cvname,
                         "beta": b, "se": se,
                         "lo": b - 1.96 * se, "hi": b + 1.96 * se,
                         "p": 2 * norm.sf(abs(b / se)) if se == se and se > 0 else np.nan,
                         "n_nights": len(d), "n_persons": d["person_id"].nunique(),
                         "n_clusters_zip3": int(cz.max() + 1),
                         "n_clusters_zip3day": int(czd.max() + 1)})
    out = pd.DataFrame(rows)
    wcb = None
    if wcb_term is not None:
        j = xs.index(wcb_term)
        se_z = np.sqrt(Vz[j, j])
        t_obs = beta[j] / se_z
        keep = [i for i in range(len(xs)) if i != j]
        Xr = X[:, keep]
        br = np.linalg.lstsq(Xr, Y, rcond=None)[0]
        ur = Y - Xr @ br
        fit0 = Xr @ br
        rng = np.random.default_rng(RNG)
        G = cz.max() + 1
        AtY = XtX_inv @ X.T          # k x n — 每次 bootstrap 只需一次 matvec
        exceed = 0
        for b_ in range(B_WCB):
            s = rng.choice([-1.0, 1.0], size=G)[cz]
            yb = fit0 + ur * s
            bb = AtY @ yb
            ub = yb - X @ bb
            Vb = cluster_cov(X, ub, cz, XtX_inv)
            tb = bb[j] / np.sqrt(Vb[j, j])
            if abs(tb) >= abs(t_obs):
                exceed += 1
        wcb = {"model": label, "term": wcb_term, "cluster": "zip3",
               "B": B_WCB, "t_obs_zip3": t_obs,
               "p_wcb": (exceed + 1) / (B_WCB + 1),
               "beta": beta[j], "se_zip3": se_z}
        log(f"  WCB {label}/{wcb_term}: t={t_obs:.3f}, p_wcb={wcb['p_wcb']:.4f}")
    return out, wcb

# ---------------------------------------------------------------- 设计变量
def heat_flag(d, q, col="tmean"):
    thr = d.groupby("zm")[col].transform(lambda s: s.quantile(q))
    return (d[col] >= thr).astype(float)

df["hd"] = heat_flag(df, 0.90)
df["sd"] = df["smoke_day"]
df["sxh"] = df["sd"] * df["hd"]
df["g_h"] = ((df.hd == 1) & (df.sd == 0)).astype(float)
df["g_s"] = ((df.hd == 0) & (df.sd == 1)).astype(float)
df["g_b"] = ((df.hd == 1) & (df.sd == 1)).astype(float)
df["tmean2"] = df["tmean"] ** 2
TS, kn_tmean = add_spline(df, "tmean", "ts")
log("tmean spline knots:", np.round(kn_tmean, 2))
NS = []
if HAVE_TMIN:
    dtm = df.dropna(subset=["tmin_night"])
    kn_tmin = np.percentile(dtm["tmin_night"], [5, 27.5, 50, 72.5, 95])
    Sm = ncs_basis(df["tmin_night"].fillna(df["tmin_night"].median()).values, kn_tmin)
    NS = [f"ns{i}" for i in range(Sm.shape[1])]
    for i, nm in enumerate(NS):
        df[nm] = np.where(df["tmin_night"].notna(), Sm[:, i], np.nan)
    log("tmin spline knots:", np.round(kn_tmin, 2))

# ================================================================ B0: 按格温度表
cell = np.select([(df.hd == 1) & (df.sd == 1), (df.hd == 1), (df.sd == 1)],
                 ["both", "hot_only", "smoke_only"], default="neither")
gb = df.assign(cell=cell).groupby("cell")
tcols = ["tmean"] + (["tmin_night"] if HAVE_TMIN else [])
parts = []
for c in tcols:
    s = gb[c]
    parts.append(pd.DataFrame({f"{c}_mean": s.mean(), f"{c}_sd": s.std(),
                               f"{c}_p50": s.quantile(0.5), f"{c}_p90": s.quantile(0.9)}))
parts.append(pd.DataFrame({"smoke_mean": gb["smoke"].mean(),
                           "tst_mean": gb["minute_asleep"].mean()}))
ct = pd.concat(parts, axis=1).join(gb.size().rename("n_nights"))
ct.loc[ct["n_nights"] < MIN_CELL] = np.nan
ct.to_csv("rev2_cell_temperature.csv")
log("=== 按格温度表 ===\n", ct.round(3).to_string())
d_tmean_both_vs_hot = (ct.loc["both", "tmean_mean"] - ct.loc["hot_only", "tmean_mean"])
log(f"both vs hot_only 的 tmean 差 = {d_tmean_both_vs_hot:.3f} C "
    f"(panel: 4C 差即可在 -0.45min/C 下伪造全部协同)")

# ================================================================ B+A: 模型阶梯
LADDER = [("M0_notemp", []),
          ("M1_tmean_lin", ["tmean"]),
          ("M2_tmean_quad", ["tmean", "tmean2"]),
          ("M3_tmean_spline", TS)]
if HAVE_TMIN:
    LADDER += [("M4_tmean_tmin_spline", TS + NS), ("M5_tmin_spline", NS)]

four_rows, int_rows, wcb_rows = [], [], []
EXCESS = {"excess_b_minus_h_minus_s": None}
for mlab, temps in LADDER:
    xs4 = ["g_h", "g_s", "g_b"] + temps
    cvec = np.zeros(len(xs4)); cvec[0], cvec[1], cvec[2] = -1, -1, 1
    o4, _ = fit_multi(df, "minute_asleep", xs4, f"four_{mlab}",
                      contrasts={"excess_b_minus_h_minus_s": cvec})
    four_rows.append(o4)
    xsi = ["sd", "hd", "sxh"] + temps
    do_wcb = mlab in ("M0_notemp", "M1_tmean_lin", "M3_tmean_spline", "M5_tmin_spline")
    oi, w = fit_multi(df, "minute_asleep", xsi, f"int_{mlab}",
                      wcb_term="sxh" if do_wcb else None)
    int_rows.append(oi)
    if w:
        wcb_rows.append(w)
    key = oi[oi["term"] == "sxh"].set_index("cov")
    log(f"{mlab}: sxh beta={key['beta'].iloc[0]:+.3f}"
        f" | p person={key.loc['person', 'p']:.4f}"
        f" zip3day={key.loc['zip3_day', 'p']:.4f}"
        f" zip3={key.loc['zip3', 'p']:.4f}"
        f" twoway={key.loc['twoway_p_z3d', 'p']:.4f}")

# 热定义换 tmin_night 的补充行 (暴露定义敏感性, 非温度调整)
if HAVE_TMIN:
    da = df.dropna(subset=["tmin_night"]).copy()
    da["hd_n"] = heat_flag(da, 0.90, col="tmin_night")
    da["sxh_n"] = da["sd"] * da["hd_n"]
    oi, w = fit_multi(da, "minute_asleep", ["sd", "hd_n", "sxh_n"] + TS,
                      "int_heatdef_tmin_night", wcb_term="sxh_n")
    int_rows.append(oi); wcb_rows.append(w)

pd.concat(four_rows).to_csv("rev2_fourcell_ladder.csv", index=False)
pd.concat(int_rows).to_csv("rev2_interaction_ladder.csv", index=False)
pd.DataFrame(wcb_rows).to_csv("rev2_wcb.csv", index=False)
log("=== 阶梯模型 + WCB 保存 ===")

# ================================================================ 3x3 重跑 (样条温度 + 多聚类)
mat = []
for sth in (10, 20, 35):
    sflag = (df["smoke"] > sth).astype(float)
    for hq, hlab in ((0.85, "P85"), (0.90, "P90"), (0.95, "P95")):
        hflag = heat_flag(df, hq)
        d = df.copy()
        d["sd_"], d["hd_"] = sflag, hflag
        d["sxh_"] = sflag * hflag
        n_b = int(d["sxh_"].sum())
        row = {"smoke_thr": sth, "heat_q": hlab, "n_both": n_b}
        if n_b < MIN_CELL:
            mat.append(row); continue
        oi, _ = fit_multi(d, "minute_asleep", ["sd_", "hd_", "sxh_"] + TS,
                          f"mx_{sth}_{hlab}")
        ii = oi[oi["term"] == "sxh_"].set_index("cov")
        row.update({"int_beta": ii.loc["person", "beta"],
                    "se_person": ii.loc["person", "se"], "p_person": ii.loc["person", "p"],
                    "se_zip3day": ii.loc["zip3_day", "se"], "p_zip3day": ii.loc["zip3_day", "p"],
                    "se_zip3": ii.loc["zip3", "se"], "p_zip3": ii.loc["zip3", "p"],
                    "se_twoway": ii.loc["twoway_p_z3d", "se"], "p_twoway": ii.loc["twoway_p_z3d", "p"]})
        mat.append(row)
pd.DataFrame(mat).to_csv("rev2_matrix3x3.csv", index=False)
log("=== 3x3 (样条温度) ===\n", pd.DataFrame(mat).round(4).to_string(index=False))

# ================================================================ C: compound x 亚组
cov = bqcsv(f"SELECT person_id, year_of_birth, has_osa, has_asthma, has_copd, "
            f"has_ihd, has_hf FROM `{DS}.person_covars`", max_rows="200000")
sub_rows = []
if cov is None and os.path.exists("person_covars_local.csv"):
    cov = pd.read_csv("person_covars_local.csv")
if cov is not None:
    d2 = df.merge(cov, on="person_id", how="inner")
    d2["cardioresp"] = (d2[["has_asthma", "has_copd", "has_ihd", "has_hf"]]
                        .sum(axis=1) > 0).astype(float)
    d2["osa_f"] = (d2["has_osa"] == 1).astype(float)
    for flag in ["cardioresp", "osa_f"]:
        for gcol in ["g_h", "g_s", "g_b"]:
            d2[f"{gcol}_x"] = d2[gcol] * d2[flag]
        xs = ["g_h", "g_s", "g_b", "g_h_x", "g_s_x", "g_b_x"] + TS
        # 亚组内 both 效应 = g_b + g_b_x; 亚组内 excess = (g_b+g_b_x)-(g_h+g_h_x)-(g_s+g_s_x)
        # excess 的组间差 = g_b_x - g_h_x - g_s_x
        k = len(xs)
        c_both_sub = np.zeros(k); c_both_sub[2] = 1; c_both_sub[5] = 1
        c_exc_sub = np.zeros(k); c_exc_sub[[2, 5]] = 1; c_exc_sub[[0, 1, 3, 4]] = -1
        c_exc_int = np.zeros(k); c_exc_int[5] = 1; c_exc_int[[3, 4]] = -1
        o, _ = fit_multi(d2, "minute_asleep", xs, f"cmpx_{flag}",
                         contrasts={"both_in_subgroup": c_both_sub,
                                    "excess_in_subgroup": c_exc_sub,
                                    "excess_interaction": c_exc_int})
        n_b1 = int(d2.loc[d2[flag] == 1, "g_b"].sum())
        n_b0 = int(d2.loc[d2[flag] == 0, "g_b"].sum())
        o["n_both_flag1"], o["n_both_flag0"] = n_b1, n_b0
        sub_rows.append(o)
        # 分层四分类 (直接可读)
        for v, tag in ((1, "yes"), (0, "no")):
            ds = d2[d2[flag] == v]
            if int(ds["g_b"].sum()) < MIN_CELL:
                continue
            cvec = np.zeros(3 + len(TS)); cvec[:3] = [-1, -1, 1]
            os_, _ = fit_multi(ds, "minute_asleep", ["g_h", "g_s", "g_b"] + TS,
                               f"strat_{flag}_{tag}",
                               contrasts={"excess_b_minus_h_minus_s": cvec})
            sub_rows.append(os_)
    pd.concat(sub_rows).to_csv("rev2_subgroup_compound.csv", index=False)
    log("=== compound x 亚组保存 ===")
else:
    log("person_covars 不可用, 跳过 C")

# ================================================================ 次级结局 (compound 格)
sec_rows = []
for y in ["minute_awake", "sleep_efficiency"]:
    if y not in df.columns or df[y].notna().sum() < 1000:
        continue
    cvec = np.zeros(3 + len(TS)); cvec[:3] = [-1, -1, 1]
    o, _ = fit_multi(df, y, ["g_h", "g_s", "g_b"] + TS, f"sec_{y}",
                     contrasts={"excess_b_minus_h_minus_s": cvec})
    sec_rows.append(o)
if sec_rows:
    pd.concat(sec_rows).to_csv("rev2_secondary_compound.csv", index=False)
    log("=== 次级结局 (compound) 保存 ===")

# ================================================================ 正主效应复查
df["tpct"] = df.groupby("zm")["tmean"].rank(pct=True)
df["smoke10c"] = df["smoke"] / 10
df["sxt"] = df["smoke10c"] * df["tpct"]
jl_rows = []
o, _ = fit_multi(df, "minute_asleep", ["smoke10c", "tpct", "sxt"], "jl_replication")
jl_rows.append(o)
o, _ = fit_multi(df, "minute_asleep", ["smoke10c", "tpct", "sxt"] + TS, "jl_tmean_spline")
jl_rows.append(o)
if HAVE_TMIN:
    o, _ = fit_multi(df, "minute_asleep", ["smoke10c", "tpct", "sxt"] + TS + NS,
                     "jl_tmean_tmin_spline")
    jl_rows.append(o)
pd.concat(jl_rows).to_csv("rev2_joint_linear_check.csv", index=False)
jj = pd.concat(jl_rows)
log("=== 连续联合模型 smoke 主效应 (person 聚类) ===\n",
    jj[(jj["term"] == "smoke10c") & (jj["cov"] == "person")]
    [["model", "beta", "lo", "hi", "p"]].round(4).to_string(index=False))

# ================================================================ 元数据
meta = pd.DataFrame([
    {"key": "rows", "value": len(df)},
    {"key": "persons", "value": df["person_id"].nunique()},
    {"key": "zip3_clusters", "value": df["zip3"].nunique()},
    {"key": "zip3day_clusters", "value": df["zip3_date"].nunique()},
    {"key": "have_tmin", "value": int(HAVE_TMIN)},
    {"key": "n_both", "value": int(df["g_b"].sum())},
    {"key": "wcb_B", "value": B_WCB},
    {"key": "tmean_knots", "value": ";".join(f"{k:.2f}" for k in kn_tmean)},
    {"key": "runtime_s", "value": int(time.time() - T0)}])
meta.to_csv("rev2_run_meta.csv", index=False)
log("ALL DONE — 导出前人工再过一遍所有 rev2_*.csv 确认无 n<20 单元格。")
log("判读顺序: (1) rev2_cell_temperature.csv 的 both-vs-hot_only tmean 差;"
    " (2) rev2_interaction_ladder.csv 里 M3/M5 的 sxh 在 zip3_day/twoway 聚类下的 p;"
    " (3) rev2_wcb.csv; (4) rev2_fourcell_ladder.csv 的 excess 行。")
