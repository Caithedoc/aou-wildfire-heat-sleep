"""19_sleep_domains_rebuild: 四域睡眠结局的数据准备 (2026-10-12, ANALYSIS_PLAN_ADDENDUM_SLEEPDOMAINS).

在 AoU Jupyter VM 里跑 (工作目录 = 放 04/12/18 脚本与 analysis_nights.parquet 的目录):
    python -u 19_sleep_domains_rebuild.py > rev8_rebuild.log 2>&1

三步, 全部 BigQuery 内聚合, 不导出个体数据:
  A. discovery: INFORMATION_SCHEMA 核 sleep_daily_summary_ext / sleep_daily_summary_counts / sleep_level
     的列名与类型 (CDR v9 支持站 2026-10 列出前两张表, RECON 0828 未见, 必须实测)。
  B. 重建 wildfire_sleep.analysis_nights: 在 04 脚本 Step 2 的 SELECT 上补齐
     minute_in_bed / minute_to_fall_asleep / minute_restless / minute_light / minute_wake,
     并挂上时相 (onset_min / offset_min / midsleep_min / midsleep_shift) 与觉醒次数 (n_awake5)。
     时相与计数优先取 _ext / _counts; 没有就从 sleep_level 聚合 (限分析样本, 谓词下推)。
  C. 删除旧 parquet 缓存并重新拉取。

列名只在 discovery 打印后按规则匹配, 匹配不到就跳过该变量并把原因写进 rev8_rebuild_meta.csv,
不猜列名 (ANALYSIS_PLAN 纪律红线 2)。
"""
import os
import re
import subprocess
import time
from io import StringIO

import pandas as pd

T0 = time.time()
CDR = "wb-silky-artichoke-2408.C2025Q4R6"          # RECON_RESULTS_20260828: Controlled Tier v9
CDRP = "wb-clever-rhubarb-4688"
DS = f"{CDRP}.wildfire_sleep"
THRESH = 20
CACHE = "analysis_nights.parquet"
meta = []


def log(*a):
    print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)


def bq(sql, quiet=False, max_rows=100000):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        f"--max_rows={max_rows}", sql], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-3000:])
    if quiet or not r.stdout.strip():
        return None
    return pd.read_csv(StringIO(r.stdout))


def note(key, value):
    meta.append({"key": key, "value": value}); log(f"meta {key} = {value}")


# ------------------------------------------------------------------ A. discovery
tabs = bq(f"""SELECT table_name FROM `{CDR}.INFORMATION_SCHEMA.TABLES`
              WHERE table_name LIKE 'sleep%' OR table_name LIKE 'heart_rate%' ORDER BY 1""")
log("sleep/heart_rate tables:\n" + tabs.to_string(index=False))
have = set(tabs["table_name"])
note("has_sleep_daily_summary_ext", "sleep_daily_summary_ext" in have)
note("has_sleep_daily_summary_counts", "sleep_daily_summary_counts" in have)
note("has_sleep_level", "sleep_level" in have)
note("hr_table", "heart_rate_intraday" if "heart_rate_intraday" in have else
     ("heart_rate_minute_level" if "heart_rate_minute_level" in have else "none"))

cols = bq(f"""SELECT table_name, column_name, data_type FROM `{CDR}.INFORMATION_SCHEMA.COLUMNS`
              WHERE table_name IN ('sleep_daily_summary','sleep_daily_summary_ext',
                                   'sleep_daily_summary_counts','sleep_level')
              ORDER BY table_name, ordinal_position""")
log("columns:\n" + cols.to_string(index=False))
cols.to_csv("rev8_sleep_table_columns.csv", index=False)


def pick(table, pattern, types=None):
    c = cols[cols.table_name == table]
    if types:
        c = c[c.data_type.isin(types)]
    hit = [x for x in c.column_name if re.search(pattern, x, re.I)]
    return hit[0] if hit else None


# ------------------------------------------------------------------ B1. 时相来源
TIMETYPES = ["TIMESTAMP", "DATETIME"]
timing_sql, timing_src = None, "none"
if "sleep_daily_summary_ext" in have:
    c_start = pick("sleep_daily_summary_ext", r"start", TIMETYPES)
    c_end = pick("sleep_daily_summary_ext", r"end", TIMETYPES)
    c_pid = pick("sleep_daily_summary_ext", r"^person_id$")
    c_date = pick("sleep_daily_summary_ext", r"sleep_date|^date$", ["DATE"])
    c_main = pick("sleep_daily_summary_ext", r"main")
    if all([c_start, c_end, c_pid, c_date]):
        main_filter = f"AND LOWER(CAST({c_main} AS STRING)) = 'true'" if c_main else ""
        timing_sql = f"""
          SELECT {c_pid} AS person_id, {c_date} AS sleep_date,
                 MIN({c_start}) AS t_onset, MAX({c_end}) AS t_offset
          FROM `{CDR}.sleep_daily_summary_ext`
          WHERE {c_pid} IN (SELECT person_id FROM `{DS}.analysis_nights`) {main_filter}
          GROUP BY 1, 2"""
        timing_src = f"sleep_daily_summary_ext ({c_start}, {c_end}, main={c_main})"
if timing_sql is None and "sleep_level" in have:
    c_pid = pick("sleep_level", r"^person_id$")
    c_date = pick("sleep_level", r"sleep_date|^date$", ["DATE"])
    c_start = pick("sleep_level", r"start", TIMETYPES)
    c_dur = pick("sleep_level", r"duration|minute")
    c_main = pick("sleep_level", r"main")
    c_level = pick("sleep_level", r"^level$|sleep_level|stage")
    if all([c_pid, c_date, c_start, c_dur]):
        main_filter = f"AND LOWER(CAST({c_main} AS STRING)) = 'true'" if c_main else ""
        timing_sql = f"""
          SELECT {c_pid} AS person_id, {c_date} AS sleep_date,
                 MIN({c_start}) AS t_onset,
                 MAX(TIMESTAMP_ADD(CAST({c_start} AS TIMESTAMP), INTERVAL CAST({c_dur} AS INT64) MINUTE)) AS t_offset
          FROM `{CDR}.sleep_level`
          WHERE {c_pid} IN (SELECT person_id FROM `{DS}.analysis_nights`) {main_filter}
          GROUP BY 1, 2"""
        timing_src = f"sleep_level ({c_start}, {c_dur}, main={c_main}, level={c_level})"
note("timing_source", timing_src)

# ------------------------------------------------------------------ B2. 觉醒次数来源
count_sql, count_src = None, "none"
if "sleep_daily_summary_counts" in have:
    c_pid = pick("sleep_daily_summary_counts", r"^person_id$")
    c_date = pick("sleep_daily_summary_counts", r"sleep_date|^date$", ["DATE"])
    c_wake = pick("sleep_daily_summary_counts", r"wake|awake")
    c_main = pick("sleep_daily_summary_counts", r"main")
    if all([c_pid, c_date, c_wake]):
        main_filter = f"AND LOWER(CAST({c_main} AS STRING)) = 'true'" if c_main else ""
        count_sql = f"""
          SELECT {c_pid} AS person_id, {c_date} AS sleep_date, MAX({c_wake}) AS n_awake_all
          FROM `{CDR}.sleep_daily_summary_counts`
          WHERE {c_pid} IN (SELECT person_id FROM `{DS}.analysis_nights`) {main_filter}
          GROUP BY 1, 2"""
        count_src = f"sleep_daily_summary_counts ({c_wake}; 注意: counts 表不区分片段长度, n_awake5 由 sleep_level 补)"
if "sleep_level" in have:
    c_pid = pick("sleep_level", r"^person_id$")
    c_date = pick("sleep_level", r"sleep_date|^date$", ["DATE"])
    c_dur = pick("sleep_level", r"duration|minute")
    c_level = pick("sleep_level", r"^level$|sleep_level|stage")
    c_main = pick("sleep_level", r"main")
    if all([c_pid, c_date, c_dur, c_level]):
        main_filter = f"AND LOWER(CAST({c_main} AS STRING)) = 'true'" if c_main else ""
        count5_sql = f"""
          SELECT {c_pid} AS person_id, {c_date} AS sleep_date,
                 COUNTIF(LOWER({c_level}) IN ('wake','awake') AND CAST({c_dur} AS FLOAT64) >= 5) AS n_awake5,
                 COUNTIF(LOWER({c_level}) IN ('wake','awake')) AS n_awake_any
          FROM `{CDR}.sleep_level`
          WHERE {c_pid} IN (SELECT person_id FROM `{DS}.analysis_nights`) {main_filter}
          GROUP BY 1, 2"""
        count_src = (count_src + " + " if count_sql else "") + f"sleep_level episodes ({c_level}, {c_dur})"
    else:
        count5_sql = None
else:
    count5_sql = None
note("awakening_source", count_src)

# ------------------------------------------------------------------ B3. 重建 analysis_nights (04 Step 2 + 补列)
bq(f"""
CREATE OR REPLACE TABLE `{DS}.analysis_nights_base` AS
WITH linked AS (
  SELECT c.person_id, c.sleep_date, c.minute_asleep, c.minute_awake, c.minute_in_bed,
         c.minute_to_fall_asleep, c.minute_restless, c.minute_light, c.minute_deep, c.minute_rem,
         c.minute_wake, c.sleep_efficiency, c.src_id,
         z.zip3,
         e.smokePM_popwt AS smoke, e.smoke_l1, e.smoke_l2, e.smoke_l3,
         e.smoke_l4, e.smoke_l5, e.smoke_l6, e.smoke_l7, e.smoke_lead1,
         e.tmean_popwt AS tmean, e.tmean_l1,
         FORMAT_DATE('%Y%m', c.sleep_date) AS ym,
         EXTRACT(DAYOFWEEK FROM c.sleep_date) AS dow
  FROM `{DS}.cohort_nights` c
  JOIN `{DS}.person_zip3` z USING (person_id)
  JOIN `{DS}.exposure_lags` e ON e.zip3 = z.zip3 AND e.date = c.sleep_date
),
strata AS (
  SELECT person_id, ym FROM linked GROUP BY person_id, ym HAVING MAX(smoke) > {THRESH}
)
SELECT l.* FROM linked l JOIN strata USING (person_id, ym)
""", quiet=True)
log("analysis_nights_base rebuilt")

parts = ["SELECT b.*"]
joins = []
if timing_sql:
    bq(f"CREATE OR REPLACE TABLE `{DS}.night_timing` AS {timing_sql}", quiet=True)
    parts.append("""
      , TIMESTAMP_DIFF(CAST(t.t_onset AS TIMESTAMP), TIMESTAMP(b.sleep_date), MINUTE) AS onset_min
      , TIMESTAMP_DIFF(CAST(t.t_offset AS TIMESTAMP), TIMESTAMP(b.sleep_date), MINUTE) AS offset_min""")
    joins.append(f"LEFT JOIN `{DS}.night_timing` t USING (person_id, sleep_date)")
    log("night_timing built")
if count_sql:
    bq(f"CREATE OR REPLACE TABLE `{DS}.night_counts` AS {count_sql}", quiet=True)
    parts.append(", k.n_awake_all")
    joins.append(f"LEFT JOIN `{DS}.night_counts` k USING (person_id, sleep_date)")
if count5_sql:
    bq(f"CREATE OR REPLACE TABLE `{DS}.night_counts5` AS {count5_sql}", quiet=True)
    parts.append(", k5.n_awake5, k5.n_awake_any")
    joins.append(f"LEFT JOIN `{DS}.night_counts5` k5 USING (person_id, sleep_date)")

bq(f"""
CREATE OR REPLACE TABLE `{DS}.analysis_nights` AS
WITH j AS (
  {' '.join(parts)}
  FROM `{DS}.analysis_nights_base` b
  {' '.join(joins)}
)
SELECT j.*
  {', (onset_min + offset_min) / 2.0 AS midsleep_min' if timing_sql else ''}
  {', (onset_min + offset_min) / 2.0 - LAG((onset_min + offset_min) / 2.0) OVER (PARTITION BY person_id ORDER BY sleep_date) AS midsleep_shift_raw' if timing_sql else ''}
  {', DATE_DIFF(sleep_date, LAG(sleep_date) OVER (PARTITION BY person_id ORDER BY sleep_date), DAY) AS gap_prev' if timing_sql else ''}
FROM j
""", quiet=True)
log("analysis_nights rebuilt with domains")

# ------------------------------------------------------------------ B4. 质量检查 (全聚合)
chk = bq(f"""SELECT COUNT(*) nights, COUNT(DISTINCT person_id) persons,
  COUNTIF(minute_deep IS NOT NULL) stage_nights,
  COUNTIF(minute_wake IS NOT NULL) stages_mode_nights,
  COUNTIF(minute_in_bed IS NOT NULL) inbed_nights,
  COUNTIF(minute_to_fall_asleep IS NOT NULL) latency_nights
  {', COUNTIF(onset_min IS NOT NULL) timing_nights' if timing_sql else ''}
  {', COUNTIF(n_awake5 IS NOT NULL) count5_nights' if count5_sql else ''}
  FROM `{DS}.analysis_nights`""")
log("QC:\n" + chk.to_string(index=False))
for k, v in chk.iloc[0].items():
    note(k, int(v))
if timing_sql:
    dist = bq(f"""SELECT CAST(FLOOR(onset_min / 60) AS INT64) AS onset_hour_bin, COUNT(*) n
                  FROM `{DS}.analysis_nights` WHERE onset_min IS NOT NULL GROUP BY 1 ORDER BY 1""")
    log("onset_min 分布 (小时桶, 相对 sleep_date 00:00; 本地时间众数应在 -3..+1):\n" + dist.to_string(index=False))
    dist.to_csv("rev8_onset_hour_dist.csv", index=False)
    mode_bin = int(dist.loc[dist.n.idxmax(), "onset_hour_bin"])
    note("onset_mode_hour_bin", mode_bin)
    note("timezone_flag", "LOCAL_OK" if -4 <= mode_bin <= 1 else "SHIFTED_CHECK_UTC")

# ------------------------------------------------------------------ C. 刷新缓存
if os.path.exists(CACHE):
    os.replace(CACHE, CACHE + ".bak_pre_rev8")
    log("old parquet moved to", CACHE + ".bak_pre_rev8")
import pandas_gbq  # noqa: E402
df = pandas_gbq.read_gbq(f"SELECT * FROM `{DS}.analysis_nights`", project_id=CDRP, progress_bar_type=None)
df.to_parquet(CACHE, index=False)
note("parquet_rows", len(df)); note("parquet_cols", ";".join(df.columns))
note("runtime_s", int(time.time() - T0))
pd.DataFrame(meta).to_csv("rev8_rebuild_meta.csv", index=False)
log("DONE. 下一步: python -u 20_sleep_domains_joint.py")
