"""22_table1_smoking: Table 1 吸烟状态行 (2026-10-30), 应 Trevathan 批注 ("Table 1 吸烟呢").

样本与四分位分组与 21_table1_quartiles.py 完全相同 (回归样本 30,268 人 / 1,633,414 夜; 人均野火 PM2.5 pd.qcut 四等分).
吸烟状态来自 All of Us Lifestyle 问卷 (PPI, CDR observation 表), 两道题:
  "100 cigarettes in lifetime" (Yes/No) 与 "smoke frequency now" (every day / some days / not at all).
  current = 100 支 Yes 且 every day 或 some days; former = Yes 且 not at all; never = No; 其余 = not reported.
每人每题取最早一次作答. 问题与答案概念按名称匹配, 匹配到的全部打印以便核对.

在 AoU Jupyter VM 跑: nohup python -u 22_table1_smoking.py > rev9_smoking.log 2>&1 &
输出 (全部聚合; 计数 1-19 置空并标 suppressed=1):
  rev9_table1_smoking.csv       group, var, cat, kind, value, value2, suppressed
  rev9_table1_smoking_meta.csv  四分位切点、人数、无问卷人数、运行信息
"""
import subprocess, time
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
CDR = "wb-silky-artichoke-2408.C2025Q4R6"
MIN_CELL = 20
T0 = time.time()


def log(*a): print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


def bqcsv(sql, max_rows="1000000"):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv", f"--max_rows={max_rows}", sql],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[:800], flush=True); return None
    return pd.read_csv(StringIO(r.stdout))


def shown(k):
    return "" if 0 < k < MIN_CELL else str(int(k))


# ---------------------------------------------------------------- 样本与四分位 (同 21)
COLS = ["person_id", "sleep_date", "ym", "dow", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["stratum"] = df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
log("rows", len(df), "| persons", df["person_id"].nunique())
assert len(df) == 1633414 and df["person_id"].nunique() == 30268, "样本与回归样本不一致"
pm = df.groupby("person_id")["smoke"].mean()
cuts = pm.quantile([0.25, 0.50, 0.75]).values
quart = pd.qcut(pm, 4, labels=["Q1", "Q2", "Q3", "Q4"]).astype(str)
log("quartile cuts", np.round(cuts, 4), "| sizes", quart.value_counts().sort_index().to_dict())

# ---------------------------------------------------------------- PPI 吸烟题目 (按名称发现)
qs = bqcsv(f"""SELECT concept_id, concept_code, concept_name, concept_class_id FROM `{CDR}.concept`
  WHERE vocabulary_id = 'PPI' AND (LOWER(concept_name) LIKE 'smoking%' OR concept_code LIKE 'Smoking_%')
  ORDER BY concept_id""", "2000")
print("\n=== PPI concepts matching smoking ===", flush=True)
print(qs.to_string(index=False), flush=True)

obs = bqcsv(f"""SELECT o.person_id, o.observation_source_concept_id AS qid, q.concept_name AS question,
       o.value_source_concept_id AS aid, a.concept_name AS answer, o.observation_date
  FROM `{CDR}.observation` o
  JOIN (SELECT DISTINCT person_id FROM `{DS}.person_covars`) s USING (person_id)
  JOIN `{CDR}.concept` q ON q.concept_id = o.observation_source_concept_id
  LEFT JOIN `{CDR}.concept` a ON a.concept_id = o.value_source_concept_id
  WHERE q.vocabulary_id = 'PPI' AND (LOWER(q.concept_name) LIKE 'smoking%' OR q.concept_code LIKE 'Smoking_%')""",
            "5000000")
assert obs is not None
obs = obs[obs["person_id"].isin(pm.index)].copy()
obs["answer"] = obs["answer"].fillna("").astype(str)
obs["observation_date"] = pd.to_datetime(obs["observation_date"], errors="coerce")
obs = obs.sort_values(["person_id", "observation_date"])
log("smoking observations among analysis persons:", len(obs), "| persons with any:", obs["person_id"].nunique())

print("\n=== question x answer counts (analysis persons; 1-19 blank) ===", flush=True)
qa = obs.groupby(["question", "answer"]).size().reset_index(name="n")
qa["n"] = qa["n"].map(shown)
print(qa.to_string(index=False), flush=True)

ql = obs["question"].str.lower()
is_hund = ql.str.contains("100")
is_freq = ql.str.contains("frequen") | ql.str.contains("every day") | ql.str.contains("now smoke")
log("question names used for 100-cig:", sorted(obs.loc[is_hund, "question"].unique()))
log("question names used for frequency:", sorted(obs.loc[is_freq, "question"].unique()))


def first_answer(mask):
    s = obs[mask].drop_duplicates("person_id", keep="first").set_index("person_id")["answer"].str.lower()
    return s.reindex(pm.index).fillna("")


hund = first_answer(is_hund)
freq = first_answer(is_freq)
ever = hund.str.contains(r"\byes\b") & ~hund.str.contains("skip|prefer|know")
never = hund.str.contains(r":\s*no\b|\bno\b") & ~hund.str.contains("skip|prefer|know|yes")
cur = freq.str.contains("every day|some days")
notatall = freq.str.contains("not at all")
status = pd.Series("Not reported", index=pm.index)
status[never] = "Never smoker"
status[ever & notatall] = "Former smoker"
status[ever & cur] = "Current smoker"
ever_unknown = int((ever & ~cur & ~notatall).sum())
no_survey = int((hund == "").sum())
log("derived status:", {k: shown(v) for k, v in status.value_counts().to_dict().items()},
    "| ever smoker with unknown current status:", shown(ever_unknown), "| no 100-cig answer at all:", shown(no_survey))

# ---------------------------------------------------------------- 按组汇总
rows = []
ORDER = ["Current smoker", "Former smoker", "Never smoker", "Not reported"]
for g in ["Overall", "Q1", "Q2", "Q3", "Q4"]:
    idx = pm.index if g == "Overall" else quart[quart == g].index
    n = len(idx)
    rows.append((g, "participants", "", "count", n, np.nan, 0))
    vc = status.reindex(idx).value_counts()
    for cat in ORDER:
        k = int(vc.get(cat, 0)); sup = int(0 < k < MIN_CELL)
        rows.append((g, "smoking", cat, "n_pct", np.nan if sup else k, np.nan if sup else k / n * 100, sup))
out = pd.DataFrame(rows, columns=["group", "var", "cat", "kind", "value", "value2", "suppressed"])
out.to_csv("rev9_table1_smoking.csv", index=False)
meta = pd.DataFrame({"key": ["pm_cut_p25", "pm_cut_p50", "pm_cut_p75", "persons", "ever_unknown_current", "no_100cig_answer",
                             "suppressed_cells", "runtime_s"],
                     "value": [cuts[0], cuts[1], cuts[2], len(pm), ever_unknown, no_survey, int(out["suppressed"].sum()),
                               round(time.time() - T0)]})
meta.to_csv("rev9_table1_smoking_meta.csv", index=False)
print("\n=== by group ===", flush=True)
print(out.to_string(index=False), flush=True)
log("done; suppressed cells:", int(out["suppressed"].sum()))
