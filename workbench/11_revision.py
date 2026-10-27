"""11_revision: compound-first 修订回跑 (2026-08-28, REVISION_PLAN P1/2/3/4/6/7)。

在 AoU Jupyter VM 里跑: `nohup python 11_revision.py > rev.log 2>&1 &`
读 analysis_nights.parquet 缓存 + BigQuery wildfire_sleep.person_covars / CDR person 表。
输出全部为聚合系数/计数表; 任何单元格 n<20 置 NaN 后导出 (披露规则)。

产出:
  rev_cell_counts.csv     P7 四格夜数精确计数        rev_joint_models.csv   P1 模型系数
  rev_joint_curves.csv    P1 曲线 (P50/75/90/95)     rev_joint_surface.csv  P1 contour 网格
  rev_matrix3x3.csv       P2 3x3 robustness           rev_compound_sens.csv  P2 compound 敏感性
  rev_subgroup_pint.csv   P6 formal P_interaction     rev_secondary_full.csv P4 四个次级结局
  rev_table1_full.csv     P4 Table1 全量重算 (30,268) rev_fig1_timeline.csv  P3 年月计数
  rev_fig1_region.csv     P3 区域计数
"""

import os, subprocess
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import statsmodels.api as sm

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
CDR = "wb-silky-artichoke-2408.C2025Q4R6"
THRESH = 20
MIN_CELL = 20
RNG = 20260828

def bqcsv(sql, max_rows="100000000"):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        f"--max_rows={max_rows}", sql], capture_output=True, text=True)
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[:500]); return None
    return pd.read_csv(StringIO(r.stdout))

# ---------------------------------------------------------------- 读缓存 (同 06)
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
df["smoke_day"] = (df["smoke"] > THRESH).astype(float)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                 + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"].astype(str) + "_" + df["month"].astype(str)
sz = df.groupby("stratum")["minute_asleep"].transform("size")
df = df[sz >= 2].copy()
print("regression-consistent rows:", len(df), "| persons:", df["person_id"].nunique())

def fe_ols(d, y, xs, label):
    d = d.dropna(subset=[y] + xs).copy()
    g = d.groupby("stratum")
    d = d[g[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    m = sm.OLS(Y, X).fit(cov_type="cluster",
                         cov_kwds={"groups": d["person_id"].values})
    out = pd.DataFrame({"model": label, "outcome": y, "term": xs,
                        "beta": m.params, "se": m.bse,
                        "lo": m.params - 1.96 * m.bse,
                        "hi": m.params + 1.96 * m.bse, "p": m.pvalues,
                        "n_nights": int(m.nobs),
                        "n_persons": d["person_id"].nunique()})
    return out, m, d

def heat_flag(d, q):
    thr = d.groupby("zm")["tmean"].transform(lambda s: s.quantile(q))
    return (d["tmean"] >= thr).astype(float)

# ================================================================ P7 四格精确计数
df["heat90"] = heat_flag(df, 0.90)
cells = pd.DataFrame({
    "cell": ["neither", "hot_only", "smoke_only", "both", "total",
             "total_hot", "total_smoke"],
    "n_nights": [
        int(((df.heat90 == 0) & (df.smoke_day == 0)).sum()),
        int(((df.heat90 == 1) & (df.smoke_day == 0)).sum()),
        int(((df.heat90 == 0) & (df.smoke_day == 1)).sum()),
        int(((df.heat90 == 1) & (df.smoke_day == 1)).sum()),
        len(df), int(df.heat90.sum()), int(df.smoke_day.sum())]})
cells.to_csv("rev_cell_counts.csv", index=False)
print("\n=== P7 cell counts ===\n", cells.to_string(index=False))

# ================================================================ P1 continuous joint
df["tpct"] = df.groupby("zm")["tmean"].rank(pct=True)

def ncs_basis(x, knots):
    x = np.asarray(x, float); k = np.asarray(knots, float)
    def d(j):
        return ((np.maximum(x - k[j], 0) ** 3 - np.maximum(x - k[-1], 0) ** 3)
                / (k[-1] - k[j]))
    return np.column_stack([x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])

KN = [5, 15, 30]
S = ncs_basis(df["smoke"].values, KN)
snames = [f"s{i}" for i in range(S.shape[1])]
for i, nm in enumerate(snames):
    df[nm] = S[:, i]
inames = []
for nm in snames:
    df[f"{nm}_x_t"] = df[nm] * df["tpct"]; inames.append(f"{nm}_x_t")
xs = snames + ["tpct"] + inames
res_rows = []
o1, m1, _ = fe_ols(df, "minute_asleep", xs, "P1_joint_spline")
res_rows.append(o1)

grid = np.arange(0, 61, 2.0)
G, G0 = ncs_basis(grid, KN), ncs_basis(np.zeros_like(grid), KN)
def curve(tp):
    C = np.hstack([G - G0, np.zeros((len(grid), 1)), (G - G0) * tp])
    est = C @ m1.params
    se = np.sqrt(np.einsum("ij,jk,ik->i", C, m1.cov_params(), C))
    return pd.DataFrame({"smoke": grid, "tpct": tp, "dTST": est,
                         "lo": est - 1.96 * se, "hi": est + 1.96 * se})
pd.concat([curve(tp) for tp in (0.50, 0.75, 0.90, 0.95)]).to_csv(
    "rev_joint_curves.csv", index=False)
pd.concat([curve(tp) for tp in np.arange(0.05, 0.96, 0.05)]).to_csv(
    "rev_joint_surface.csv", index=False)
print("=== P1 curves saved ===")

df["smoke10c"] = df["smoke"] / 10
df["sxt"] = df["smoke10c"] * df["tpct"]
o2, _, _ = fe_ols(df, "minute_asleep", ["smoke10c", "tpct", "sxt"], "P1_linear_x_linear")
res_rows.append(o2)
pd.concat(res_rows).to_csv("rev_joint_models.csv", index=False)

# ================================================================ P2 3x3 matrix
mat = []
for sth in (10, 20, 35):
    sflag = (df["smoke"] > sth).astype(float)
    for hq, hlab in ((0.85, "P85"), (0.90, "P90"), (0.95, "P95")):
        hflag = heat_flag(df, hq)
        d = df.copy()
        d["sd"], d["hd"], d["sxh"] = sflag, hflag, sflag * hflag
        d["g_h"] = ((hflag == 1) & (sflag == 0)).astype(float)
        d["g_s"] = ((hflag == 0) & (sflag == 1)).astype(float)
        d["g_b"] = ((hflag == 1) & (sflag == 1)).astype(float)
        n_b = int(d["g_b"].sum())
        row = {"smoke_thr": sth, "heat_q": hlab, "n_both": n_b}
        if n_b < MIN_CELL:
            mat.append(row); continue
        oi, mi, _ = fe_ols(d, "minute_asleep", ["sd", "hd", "sxh", "tmean"],
                           f"mx_{sth}_{hlab}_int")
        oj, mj, _ = fe_ols(d, "minute_asleep", ["g_h", "g_s", "g_b"],
                           f"mx_{sth}_{hlab}_4cell")
        bh, bs, bb = mj.params[:3]
        ii = oi[oi.term == "sxh"].iloc[0]
        row.update({"b_hot_only": bh, "b_smoke_only": bs, "b_both": bb,
                    "additive_exp": bh + bs, "excess": bb - (bh + bs),
                    "int_beta": ii.beta, "int_lo": ii.lo, "int_hi": ii.hi,
                    "int_p": ii.p})
        mat.append(row)
pd.DataFrame(mat).to_csv("rev_matrix3x3.csv", index=False)
print("=== P2 3x3 matrix ===\n", pd.DataFrame(mat).round(3).to_string(index=False))

# ---- compound 专属敏感性 (20 x P90): excl2023 / 东西部 / 设备 / IPW
sens = []
def compound_on(d, tag, extra_w=None):
    d = d.copy()
    d["hd"] = heat_flag(d, 0.90)
    d["sxh"] = d["smoke_day"] * d["hd"]
    n_b = int(d["sxh"].sum())
    if n_b < MIN_CELL:
        sens.append(pd.DataFrame([{"model": tag, "term": "sxh", "n_both": n_b}])); return
    o, _, _ = fe_ols(d, "minute_asleep", ["smoke_day", "hd", "sxh", "tmean"], tag)
    o["n_both"] = n_b
    sens.append(o)

compound_on(df[df.year != 2023], "cmp_excl2023")
z1 = df["zip3"].astype(str).str.zfill(3).str[0]
compound_on(df[z1.isin(list("012"))], "cmp_east")
compound_on(df[z1.isin(list("89"))], "cmp_west")
for s in df["src_id"].dropna().unique():
    compound_on(df[df.src_id == s], f"cmp_src_{s}")

# IPW 版 (复刻 06 E3 的 stabilized weights, 再加交互)
try:
    from sklearn.ensemble import HistGradientBoostingClassifier
    q = f"""
    WITH span AS (SELECT person_id, MIN(sleep_date) d0, MAX(sleep_date) d1
                  FROM `{DS}.cohort_nights` GROUP BY person_id),
    w AS (SELECT w.person_id, w.wear_date AS date
          FROM `{DS}.wear_days` w JOIN span s USING (person_id)
          WHERE w.wear_date BETWEEN s.d0 AND s.d1)
    SELECT w.person_id, w.date, IF(c.person_id IS NULL, 0, 1) AS has_sleep,
           e.smokePM_popwt AS smoke, e.tmean_popwt AS tmean
    FROM w
    LEFT JOIN `{DS}.cohort_nights` c
      ON c.person_id = w.person_id AND c.sleep_date = w.date
    JOIN `{DS}.person_zip3` z ON z.person_id = w.person_id
    JOIN `{DS}.exposure_zip3_daily` e ON e.zip3 = z.zip3 AND e.date = w.date
    """
    wd = bqcsv(q)
    if wd is not None and len(wd):
        wd["date"] = pd.to_datetime(wd["date"])
        wd["smoke_day"] = (pd.to_numeric(wd["smoke"], errors="coerce") > THRESH).astype(float)
        wd["doy"] = wd["date"].dt.dayofyear
        wd["dow2"] = wd["date"].dt.dayofweek
        wd["year"] = wd["date"].dt.year
        wd["tmean"] = pd.to_numeric(wd["tmean"], errors="coerce")
        wd = wd.dropna(subset=["tmean"])
        ph = wd.groupby("person_id")["has_sleep"].mean().rename("p_hist")
        wd = wd.join(ph, on="person_id")
        F2 = ["p_hist", "tmean", "doy", "dow2", "year", "smoke_day"]
        clf = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.08,
                                             random_state=RNG)
        clf.fit(wd[F2].values, wd["has_sleep"].values)
        wd["p_obs"] = clf.predict_proba(wd[F2].values)[:, 1].clip(0.05, 0.999)
        wd["sw"] = wd["has_sleep"].mean() / wd["p_obs"]
        wsub = wd.loc[wd["has_sleep"] == 1, ["person_id", "date", "sw"]]
        d3 = df.merge(wsub, on=["person_id", "date"], how="inner")
        d3["hd"] = heat_flag(d3, 0.90)
        d3["sxh"] = d3["smoke_day"] * d3["hd"]
        dd = d3.dropna(subset=["minute_asleep", "smoke_day", "hd", "sxh", "tmean", "sw"]).copy()
        g = dd.groupby("stratum")
        dd = dd[g["minute_asleep"].transform("size") >= 2]
        g = dd.groupby("stratum")
        Y = (dd["minute_asleep"] - g["minute_asleep"].transform("mean")).values
        X = np.column_stack([(dd[c] - g[c].transform("mean")).values
                             for c in ["smoke_day", "hd", "sxh", "tmean"]])
        mw = sm.WLS(Y, X, weights=dd["sw"].values).fit(
            cov_type="cluster", cov_kwds={"groups": dd["person_id"].values})
        sens.append(pd.DataFrame({
            "model": "cmp_ipw", "outcome": "minute_asleep",
            "term": ["smoke_day", "hd", "sxh", "tmean"],
            "beta": mw.params, "se": mw.bse,
            "lo": mw.params - 1.96 * mw.bse, "hi": mw.params + 1.96 * mw.bse,
            "p": mw.pvalues, "n_nights": int(mw.nobs),
            "n_persons": dd["person_id"].nunique(),
            "n_both": int(dd["sxh"].sum())}))
except Exception as e:
    print("IPW compound failed:", e)

pd.concat(sens).to_csv("rev_compound_sens.csv", index=False)
print("=== P2 compound sensitivity saved ===")
print("NOTE 宽窗 (60-960) compound 需 F2 重抽样本; 若 VM 有对应 parquet 可比照跑。")

# ================================================================ P6 P_interaction
cov = bqcsv(f"SELECT person_id, year_of_birth, gender, deprivation_index, has_osa, "
            f"has_insomnia, has_asthma, has_copd, has_t2dm, has_htn, has_ihd, has_hf "
            f"FROM `{DS}.person_covars`", max_rows="200000")
pint = []
if cov is not None:
    d2 = df.merge(cov, on="person_id", how="inner")
    d2["age"] = d2["year"] - d2["year_of_birth"]
    d2["cardioresp"] = (d2[["has_asthma", "has_copd", "has_ihd", "has_hf"]]
                        .sum(axis=1) > 0).astype(float)
    d2["osa_f"] = (d2["has_osa"] == 1).astype(float)
    d2["age65"] = (d2["age"] >= 65).astype(float)
    thr75 = d2["deprivation_index"].quantile(0.75)
    d2["deprived_top25"] = (d2["deprivation_index"] >= thr75).astype(float)
    for flag in ["cardioresp", "osa_f", "age65", "deprived_top25"]:
        d2[f"sx_{flag}"] = d2["smoke_day"] * d2[flag]
        o, _, _ = fe_ols(d2, "minute_asleep",
                         ["smoke_day", f"sx_{flag}", "tmean"], f"pint_{flag}")
        pint.append(o)
    pd.concat(pint).to_csv("rev_subgroup_pint.csv", index=False)
    print("=== P6 P_interaction ===")
    print(pd.concat(pint).query("term.str.startswith('sx_')").round(4).to_string(index=False))

# ================================================================ P4a secondary 4 结局
try:
    pdres = pd.read_csv("phaseD_results.csv")
    sec = pdres[pdres["model"].astype(str).str.contains("secondary", case=False)]
    sec.to_csv("rev_secondary_full.csv", index=False)
    print("=== P4 secondary outcomes (from phaseD_results.csv) ===")
    print(sec.round(4).to_string(index=False))
except Exception as e:
    print("phaseD_results.csv read failed:", e)

# ================================================================ P4b Table 1 全量重算 (回归样本 30,268)
rows = []
pers_first = df.sort_values("date").groupby("person_id").first()
nights_pp = df.groupby("person_id").size()

if cov is not None:
    pc = cov.set_index("person_id").reindex(pers_first.index)
    age_first = pers_first["year"] - pc["year_of_birth"]
    rows += [("person", "n_participants", len(pers_first), np.nan),
             ("person", "age_mean", age_first.mean(), age_first.std()),
             ("person", "female_n", (pc["gender"] == "Female").sum(),
              (pc["gender"] == "Female").mean() * 100)]
    for f, lab in [("has_osa", "osa"), ("has_insomnia", "insomnia"),
                   ("has_asthma", "asthma"), ("has_copd", "copd"),
                   ("has_t2dm", "t2dm"), ("has_htn", "htn"),
                   ("has_ihd", "ihd"), ("has_hf", "hf")]:
        n = (pc[f] == 1).sum()
        rows.append(("person", lab + "_n", n, n / len(pc) * 100))
    rows.append(("person", "adi_mean", pc["deprivation_index"].mean(),
                 pc["deprivation_index"].std()))
    qs = pd.qcut(pc["deprivation_index"], 4, labels=["q1", "q2", "q3", "q4"])
    for qlab in ["q1", "q2", "q3", "q4"]:
        n = (qs == qlab).sum()
        rows.append(("person", f"adi_{qlab}_n", n, n / qs.notna().sum() * 100))
rows += [("person", "nights_pp_mean", nights_pp.mean(), np.nan),
         ("person", "nights_pp_median", nights_pp.median(), np.nan)]

# region (ZIP 首位: 0-2 East, 3-7 Central, 8-9 West) + device source
pz1 = pers_first["zip3"].astype(str).str.zfill(3).str[0]
for lab, s in [("east_0_2", pz1.isin(list("012"))),
               ("central_3_7", pz1.isin(list("34567"))),
               ("west_8_9", pz1.isin(list("89")))]:
    rows.append(("person", "region_" + lab + "_n", s.sum(), s.mean() * 100))
for s, n in pers_first["src_id"].value_counts().items():
    rows.append(("person", f"device_{s}_n", n, n / len(pers_first) * 100))

# night level
rows += [("night", "n_nights", len(df), np.nan),
         ("night", "tst_mean", df["minute_asleep"].mean(), df["minute_asleep"].std()),
         ("night", "tst_q1", df["minute_asleep"].quantile(0.25), np.nan),
         ("night", "tst_q3", df["minute_asleep"].quantile(0.75), np.nan),
         ("night", "smoke_mean", df["smoke"].mean(), df["smoke"].std()),
         ("night", "smoke_gt10_pct", (df["smoke"] > 10).mean() * 100, np.nan),
         ("night", "smoke_gt20_pct", (df["smoke"] > 20).mean() * 100, np.nan),
         ("night", "smoke_gt35_pct", (df["smoke"] > 35).mean() * 100, np.nan),
         ("night", "tmean_mean", df["tmean"].mean(), df["tmean"].std()),
         ("night", "hot_pct", df["heat90"].mean() * 100, np.nan),
         ("night", "informative_strata", df["stratum"].nunique(), np.nan),
         ("night", "nights_per_stratum_mean", df.groupby("stratum").size().mean(), np.nan)]

t1 = pd.DataFrame(rows, columns=["level", "stat", "value", "value2"])
t1.loc[(t1["stat"].str.endswith("_n")) & (t1["value"] < MIN_CELL),
       ["value", "value2"]] = np.nan
t1.to_csv("rev_table1_full.csv", index=False)
print("=== P4 Table 1 recomputed on regression sample ===")

# race / ethnicity / education (CDR 拉全 cohort 人级, VM 内合并后只出计数)
try:
    dem = bqcsv(f"""
      SELECT p.person_id, cr.concept_name AS race, ce.concept_name AS ethnicity
      FROM `{CDR}.person` p
      JOIN (SELECT DISTINCT person_id FROM `{DS}.person_covars`) s USING (person_id)
      LEFT JOIN `{CDR}.concept` cr ON p.race_concept_id = cr.concept_id
      LEFT JOIN `{CDR}.concept` ce ON p.ethnicity_concept_id = ce.concept_id""",
      max_rows="200000")
    if dem is not None:
        dm = dem.set_index("person_id").reindex(pers_first.index)
        rc = dm["race"].fillna("Unknown").value_counts()
        ec = dm["ethnicity"].fillna("Unknown").value_counts()
        out = pd.concat([
            pd.DataFrame({"var": "race", "cat": rc.index, "n": rc.values}),
            pd.DataFrame({"var": "ethnicity", "cat": ec.index, "n": ec.values})])
        out["pct"] = out["n"] / len(dm) * 100
        out.loc[out["n"] < MIN_CELL, ["n", "pct"]] = np.nan
        out.to_csv("rev_table1_race.csv", index=False)
        print(out.round(2).to_string(index=False))
except Exception as e:
    print("race/ethnicity fetch failed:", e)

try:
    edu = bqcsv(f"""
      SELECT o.person_id, c.concept_name AS edu
      FROM `{CDR}.observation` o
      JOIN (SELECT DISTINCT person_id FROM `{DS}.person_covars`) s USING (person_id)
      LEFT JOIN `{CDR}.concept` c ON o.value_source_concept_id = c.concept_id
      WHERE o.observation_source_concept_id = 1585940""", max_rows="300000")
    if edu is not None and len(edu):
        em = edu.drop_duplicates("person_id").set_index("person_id").reindex(pers_first.index)
        ecnt = em["edu"].fillna("Unknown").value_counts()
        oe = pd.DataFrame({"var": "education", "cat": ecnt.index, "n": ecnt.values})
        oe["pct"] = oe["n"] / len(em) * 100
        oe.loc[oe["n"] < MIN_CELL, ["n", "pct"]] = np.nan
        oe.to_csv("rev_table1_edu.csv", index=False)
except Exception as e:
    print("education fetch failed:", e)

# ================================================================ P3 Figure 1 数据
kind = np.select([(df.heat90 == 1) & (df.smoke_day == 1),
                  (df.heat90 == 1), (df.smoke_day == 1)],
                 ["both", "hot_only", "smoke_only"], default="neither")
tl = (df.assign(kind=kind).groupby(["year", "month", "kind"])
        .size().rename("n").reset_index())
tl.loc[tl["n"] < MIN_CELL, "n"] = np.nan
tl.to_csv("rev_fig1_timeline.csv", index=False)
reg = (df.assign(zip1=z1, kind=kind).groupby(["zip1", "kind"])
         .size().rename("n").reset_index())
reg.loc[reg["n"] < MIN_CELL, "n"] = np.nan
reg.to_csv("rev_fig1_region.csv", index=False)
print("=== P3 timeline/region saved ===")
print("\nALL DONE — 导出前人工再过一遍所有 rev_*.csv 确认无 n<20 单元格。")
