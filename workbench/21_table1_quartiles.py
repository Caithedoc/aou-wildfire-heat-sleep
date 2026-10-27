"""21_table1_quartiles: Table 1 按参与者平均野火 PM2.5 四分位分列 (2026-10-13).

仿 Luo 2025 JAMA Netw Open (锂暴露五分位) 与 Luo 2024 (Total + 分层列): 列 = Overall + Q1-Q4.
分层变量 = 每位参与者在其全部分析夜上的平均 wildfire-specific PM2.5 (人级, pd.qcut 四等分).
样本与 15/20 脚本完全相同 (回归样本 30,268 人 / 1,633,414 夜); Overall 列须复现现行 Table 1.

在 AoU Jupyter VM 跑: nohup python -u 21_table1_quartiles.py > rev9.log 2>&1 &
需要: analysis_nights.parquet; bq 可用 (person_covars + CDR person/observation).
输出 (全部聚合; 计数 1-20 置空并标 suppressed=1):
  rev9_table1_quartiles.csv  group, section, var, cat, kind, value, value2, suppressed
  rev9_table1_meta.csv       四分位切点 (人均 PM2.5 的 P25/P50/P75), 每组人数/夜数, 运行信息
"""
import os, subprocess, time
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
CDR = "wb-silky-artichoke-2408.C2025Q4R6"
THRESH, MIN_CELL = 20, 20
T0 = time.time()


def log(*a): print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


def bqcsv(sql, max_rows="1000000"):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv", f"--max_rows={max_rows}", sql],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[:500]); return None
    return pd.read_csv(StringIO(r.stdout))


# ---------------------------------------------------------------- 样本 (同 15/20)
COLS = ["person_id", "sleep_date", "ym", "dow", "src_id", "zip3", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["year"] = df["date"].dt.year
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["date"].dt.month.astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
log("rows", len(df), "| persons", df["person_id"].nunique())
assert len(df) == 1633414 and df["person_id"].nunique() == 30268, "样本与回归样本不一致"

df["sd"] = (df["smoke"] > THRESH).astype(int)
df["hd"] = (df["tmean"] >= df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(int)
df["both"] = df["sd"] * df["hd"]

# ---------------------------------------------------------------- 人级分层: 平均野火 PM2.5 四分位
pm = df.groupby("person_id")["smoke"].mean()
cuts = pm.quantile([0.25, 0.50, 0.75]).values
quart = pd.qcut(pm, 4, labels=["Q1", "Q2", "Q3", "Q4"])
df["grp"] = df["person_id"].map(quart).astype(str)
first = df.sort_values("date").groupby("person_id").first()
first["grp"] = quart.reindex(first.index).astype(str)

# ---------------------------------------------------------------- 人级协变量
cov = bqcsv(f"SELECT person_id, year_of_birth, gender, deprivation_index, has_osa, has_insomnia, has_asthma, "
            f"has_copd, has_t2dm, has_htn, has_ihd, has_hf FROM `{DS}.person_covars`", max_rows="200000")
assert cov is not None
pc = cov.drop_duplicates("person_id").set_index("person_id").reindex(first.index)
dem = bqcsv(f"""SELECT p.person_id, cr.concept_name AS race, ce.concept_name AS ethnicity
  FROM `{CDR}.person` p JOIN (SELECT DISTINCT person_id FROM `{DS}.person_covars`) s USING (person_id)
  LEFT JOIN `{CDR}.concept` cr ON p.race_concept_id = cr.concept_id
  LEFT JOIN `{CDR}.concept` ce ON p.ethnicity_concept_id = ce.concept_id""", max_rows="200000")
edu = bqcsv(f"""SELECT o.person_id, c.concept_name AS edu FROM `{CDR}.observation` o
  JOIN (SELECT DISTINCT person_id FROM `{DS}.person_covars`) s USING (person_id)
  LEFT JOIN `{CDR}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id = 1585940""", max_rows="300000")
dm = dem.drop_duplicates("person_id").set_index("person_id").reindex(first.index)
em = edu.drop_duplicates("person_id").set_index("person_id").reindex(first.index)


def race_cat(x):
    x = "" if pd.isna(x) else str(x)
    if x == "White": return "White"
    if x == "Black or African American": return "Black or African American"
    if x in ("Asian", "Native Hawaiian or Other Pacific Islander"): return "Asian, Native Hawaiian, or Other Pacific Islander"
    if x == "More than one population": return "More than one population"
    if x in ("Middle Eastern or North African", "None of these"): return "Middle Eastern or North African, or none of these"
    return "Not indicated or skipped"


def edu_cat(x):
    x = "" if pd.isna(x) else str(x)
    if "Advanced Degree" in x: return "Advanced degree"
    if "College Graduate" in x: return "College graduate"
    if "College One to Three" in x: return "Some college"
    if "Twelve Or GED" in x: return "High school diploma or GED"
    return "Less than high school or not reported"


P = pd.DataFrame(index=first.index)
P["grp"] = first["grp"]
P["age"] = first["year"] - pc["year_of_birth"]
P["female"] = (pc["gender"] == "Female").astype(int)
P["race"] = dm["race"].map(race_cat)
P["hisp"] = (dm["ethnicity"] == "Hispanic or Latino").astype(int)
P["edu"] = em["edu"].map(edu_cat)
z1 = first["zip3"].str[0]
P["region"] = np.select([z1.isin(list("012")), z1.isin(list("34567")), z1.isin(list("89"))],
                        ["East (0xx-2xx)", "Central/South (3xx-7xx)", "West (8xx-9xx)"], "Unknown")
P["device"] = first["src_id"].astype(str)
for f in ["has_osa", "has_insomnia", "has_asthma", "has_copd", "has_t2dm", "has_htn", "has_ihd", "has_hf"]:
    P[f] = (pc[f] == 1).astype(int)
P["adi"] = pc["deprivation_index"]
P["nights"] = df.groupby("person_id").size().reindex(P.index)
P["pm_mean"] = pm.reindex(P.index)

# ---------------------------------------------------------------- 汇总
rows = []


def add(g, section, var, cat, kind, v, v2=np.nan):
    sup = int(kind == "n_pct" and 0 < v < MIN_CELL)
    if sup: v, v2 = np.nan, np.nan
    rows.append((g, section, var, cat, kind, v, v2, sup))


ORDER = {"race": ["White", "Black or African American", "Asian, Native Hawaiian, or Other Pacific Islander", "More than one population",
                  "Middle Eastern or North African, or none of these", "Not indicated or skipped"],
         "edu": ["Advanced degree", "College graduate", "Some college", "High school diploma or GED", "Less than high school or not reported"],
         "region": ["East (0xx-2xx)", "Central/South (3xx-7xx)", "West (8xx-9xx)"]}
for g in ["Overall", "Q1", "Q2", "Q3", "Q4"]:
    pp = P if g == "Overall" else P[P["grp"] == g]
    nn = df if g == "Overall" else df[df["grp"] == g]
    n = len(pp)
    add(g, "participant", "participants", "", "count", n)
    add(g, "participant", "age", "", "mean_sd", pp["age"].mean(), pp["age"].std())
    add(g, "participant", "female", "", "n_pct", pp["female"].sum(), pp["female"].mean() * 100)
    for var in ["race", "edu", "region"]:
        vc = pp[var].value_counts()
        for cat in ORDER[var]:
            k = int(vc.get(cat, 0)); add(g, "participant", var, cat, "n_pct", k, k / n * 100)
        other = int(n - sum(vc.get(c, 0) for c in ORDER[var]))
        if other: add(g, "participant", var, "Unknown", "n_pct", other, other / n * 100)
    add(g, "participant", "hispanic", "", "n_pct", pp["hisp"].sum(), pp["hisp"].mean() * 100)
    for dev, k in pp["device"].value_counts().items():
        add(g, "participant", "device", dev, "n_pct", int(k), k / n * 100)
    for f in ["has_osa", "has_insomnia", "has_asthma", "has_copd", "has_t2dm", "has_htn", "has_ihd", "has_hf"]:
        add(g, "participant", f, "", "n_pct", pp[f].sum(), pp[f].mean() * 100)
    add(g, "participant", "adi", "", "mean_sd", pp["adi"].mean(), pp["adi"].std())
    add(g, "participant", "nights_pp", "median", "median_iqr", pp["nights"].median(), np.nan)
    add(g, "participant", "nights_pp", "q1q3", "q1q3", pp["nights"].quantile(0.25), pp["nights"].quantile(0.75))
    add(g, "participant", "pm_person_mean", "", "median_iqr", pp["pm_mean"].median(), np.nan)
    add(g, "participant", "pm_person_mean", "q1q3", "q1q3", pp["pm_mean"].quantile(0.25), pp["pm_mean"].quantile(0.75))
    m = len(nn)
    add(g, "night", "nights", "", "count", m)
    add(g, "night", "tst", "", "mean_sd", nn["minute_asleep"].mean(), nn["minute_asleep"].std())
    add(g, "night", "pm", "", "mean_sd", nn["smoke"].mean(), nn["smoke"].std())
    add(g, "night", "smoke_night", "", "n_pct", nn["sd"].sum(), nn["sd"].mean() * 100)
    add(g, "night", "hot_night", "", "n_pct", nn["hd"].sum(), nn["hd"].mean() * 100)
    add(g, "night", "both_night", "", "n_pct", nn["both"].sum(), nn["both"].mean() * 100)
    add(g, "night", "tmean", "", "mean_sd", nn["tmean"].mean(), nn["tmean"].std())

out = pd.DataFrame(rows, columns=["group", "section", "var", "cat", "kind", "value", "value2", "suppressed"])
out.to_csv("rev9_table1_quartiles.csv", index=False)
meta = pd.DataFrame({"key": ["pm_cut_p25", "pm_cut_p50", "pm_cut_p75", "persons", "nights", "suppressed_cells", "runtime_s"],
                     "value": [cuts[0], cuts[1], cuts[2], len(P), len(df), int(out["suppressed"].sum()), round(time.time() - T0)]})
meta.to_csv("rev9_table1_meta.csv", index=False)
log("done; suppressed cells:", int(out["suppressed"].sum()))
print(out[out.group == "Overall"].to_string(index=False))
