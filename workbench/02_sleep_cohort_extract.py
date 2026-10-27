"""Phase B: 队列与结局提取。在 aou-wildfire-sleep workspace 的 AoU Jupyter 里跑。

前提 (2026-08-28 Phase A 实测): Controlled Tier CDR = wb-silky-artichoke-2408.C2025Q4R6;
表名列名全部经 INFORMATION_SCHEMA 核实, 见 RECON_RESULTS_20260828.md。

产出: workspace 项目下建 dataset `wildfire_sleep`, 落四张中间表。
只落参与者级中间表在仓内, 不导出任何行级数据。

用法: 逐 cell 贴进 notebook。先跑 Cell 0 建 dataset。
"""

# %% Cell 0: 环境
import subprocess

import pandas as pd

CDR = "wb-silky-artichoke-2408.C2025Q4R6"
PROJ = "wb-clever-rhubarb-4688"          # workspace GCP project (实测)
OUT = f"{PROJ}.wildfire_sleep"           # 输出 dataset
WIN_START, WIN_END = "2016-01-01", "2023-12-31"  # 主分析窗 (暴露端封顶)


def bq(sql: str):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv", sql],
                       capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr)
    from io import StringIO
    return pd.read_csv(StringIO(r.stdout)) if r.stdout.strip() else None


# pet SA 无 bigquery.datasets.create 权限 (2026-08-28 实测), 必须让 Workbench 代建:
#   wb resource create bq-dataset --id=wildfire_sleep
# 建好后 OUT 即 wb-clever-rhubarb-4688.wildfire_sleep。
subprocess.run(["wb", "resource", "create", "bq-dataset", "--id=wildfire_sleep"],
               capture_output=True)
print("dataset ready:", OUT)

# %% Cell 1: 逐夜结局表 cohort_nights
# 夜级清理规则 (ANALYSIS_PLAN 数据清理规范): 主睡眠段; 3-14h 可信窗;
# 负值/效率>100 剔除; 同人同日取 minute_asleep 最长一条; 起始日挂接。
bq(f"""
CREATE OR REPLACE TABLE `{OUT}.cohort_nights` AS
WITH main_sleep AS (
  SELECT person_id, sleep_date, minute_in_bed, minute_asleep, minute_awake,
         minute_restless, minute_deep, minute_light, minute_rem, minute_wake,
         minute_to_fall_asleep, src_id,
         ROW_NUMBER() OVER (PARTITION BY person_id, sleep_date
                            ORDER BY minute_asleep DESC) AS rk
  FROM `{CDR}.sleep_daily_summary`
  WHERE LOWER(is_main_sleep) = 'true'
    AND sleep_date BETWEEN '{WIN_START}' AND '{WIN_END}'
)
SELECT person_id, sleep_date,
       minute_asleep, minute_in_bed, minute_awake, minute_restless,
       minute_deep, minute_light, minute_rem, minute_wake, minute_to_fall_asleep,
       SAFE_DIVIDE(minute_asleep, NULLIF(minute_in_bed, 0)) AS sleep_efficiency,
       src_id
FROM main_sleep
WHERE rk = 1
  AND minute_asleep BETWEEN 180 AND 840          -- 3-14 h 可信窗 (敏感性 60-960)
  AND minute_in_bed >= minute_asleep             -- 效率 <= 100%
  AND minute_in_bed > 0
  AND minute_asleep >= 0 AND minute_awake >= 0
""")
print(bq(f"SELECT COUNT(*) n, COUNT(DISTINCT person_id) p FROM `{OUT}.cohort_nights`"))

# %% Cell 2: person 级 zip3 表 person_zip3 (含搬家标记)
# zip 概念 1585250 实测覆盖 62,906/62,911。Controlled Tier 的 value_as_string
# 格式先抽样看一眼再定截取规则 (常见为 3 位 zip3 或全 zip; 此处按前 3 位取)。
print(bq(f"""
SELECT LENGTH(TRIM(value_as_string)) AS len, COUNT(*) AS n
FROM `{CDR}.observation`
WHERE observation_source_concept_id = 1585250
GROUP BY len ORDER BY n DESC LIMIT 5
"""))

# 2026-08-28 实测: zip 值全部为 "###**" 掩码 (746,980 条, 纯数字 0 条), 且每人只有
# 一条快照记录, 所以 moved_flag 恒为 False。CDR 无地址史, 搬家敏感性改写为局限性;
# case-crossover 同月自身对照天然压制搬家错分。
# 注意 BigQuery 分析函数不支持 COUNT(DISTINCT ...) OVER, 用 GROUP BY 改写。
bq(f"""
CREATE OR REPLACE TABLE `{OUT}.person_zip3` AS
WITH zz AS (
  SELECT person_id, SUBSTR(TRIM(value_as_string), 1, 3) AS zip3, observation_date
  FROM `{CDR}.observation`
  WHERE observation_source_concept_id = 1585250 AND value_as_string IS NOT NULL
),
latest AS (
  SELECT person_id, zip3,
         ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY observation_date DESC) AS rk
  FROM zz
),
nz AS (SELECT person_id, COUNT(DISTINCT zip3) AS n FROM zz GROUP BY person_id)
SELECT l.person_id, l.zip3, nz.n > 1 AS moved_flag
FROM latest l JOIN nz USING (person_id)
WHERE l.rk = 1
""")
print(bq(f"SELECT COUNT(*) n, COUNTIF(moved_flag) movers FROM `{OUT}.person_zip3`"))
# as-run 2026-08-28: n=746,980, movers=0

# %% Cell 3: person 级协变量表 person_covars
# 人口学 + zip3 SES + EHR 基线慢病旗标 (SNOMED 标准概念 + 后代)。
# 慢病概念根: OSA 313459, 失眠 436962, 哮喘 317009, COPD 255573,
# 2 型糖尿病 201826, 高血压 316866, 缺血性心脏病 4185932, 心衰 316139。
# (标准 OMOP 根概念; 跑前用 concept 表抽查名称对得上再继续。)
print(bq(f"""
SELECT concept_id, concept_name FROM `{CDR}.concept`
WHERE concept_id IN (313459, 436962, 317009, 255573, 201826, 316866, 4185932, 316139)
"""))

bq(f"""
CREATE OR REPLACE TABLE `{OUT}.person_covars` AS
WITH conds AS (
  SELECT co.person_id, ca.ancestor_concept_id AS root_id
  FROM `{CDR}.condition_occurrence` co
  JOIN `{CDR}.concept_ancestor` ca
    ON co.condition_concept_id = ca.descendant_concept_id
  WHERE ca.ancestor_concept_id IN (313459, 436962, 317009, 255573,
                                   201826, 316866, 4185932, 316139)
  GROUP BY co.person_id, ca.ancestor_concept_id
)
SELECT p.person_id,
       p.year_of_birth,
       c_gender.concept_name  AS gender,
       z.zip3, z.moved_flag,
       ses.deprivation_index, ses.median_income, ses.fraction_poverty,
       MAX(IF(conds.root_id = 313459, 1, 0))  AS has_osa,
       MAX(IF(conds.root_id = 436962, 1, 0))  AS has_insomnia,
       MAX(IF(conds.root_id = 317009, 1, 0))  AS has_asthma,
       MAX(IF(conds.root_id = 255573, 1, 0))  AS has_copd,
       MAX(IF(conds.root_id = 201826, 1, 0))  AS has_t2dm,
       MAX(IF(conds.root_id = 316866, 1, 0))  AS has_htn,
       MAX(IF(conds.root_id = 4185932, 1, 0)) AS has_ihd,
       MAX(IF(conds.root_id = 316139, 1, 0))  AS has_hf
FROM `{CDR}.person` p
JOIN (SELECT DISTINCT person_id FROM `{OUT}.cohort_nights`) s USING (person_id)
LEFT JOIN `{CDR}.concept` c_gender ON p.gender_concept_id = c_gender.concept_id
LEFT JOIN `{OUT}.person_zip3` z USING (person_id)
LEFT JOIN `{CDR}.zip3_ses_map` ses ON z.zip3 = ses.zip3_as_string
LEFT JOIN conds ON p.person_id = conds.person_id
GROUP BY p.person_id, p.year_of_birth, gender, z.zip3, z.moved_flag,
         ses.deprivation_index, ses.median_income, ses.fraction_poverty
""")
print(bq(f"SELECT COUNT(*) n FROM `{OUT}.person_covars`"))

# %% Cell 4: 佩戴日历表 wear_days (IPW 打底)
# 佩戴日 = 当日 heart_rate_summary 或 steps_intraday 有任意记录。
bq(f"""
CREATE OR REPLACE TABLE `{OUT}.wear_days` AS
SELECT person_id, dt AS wear_date
FROM (
  SELECT DISTINCT person_id, date AS dt
  FROM `{CDR}.heart_rate_summary`                 -- 日级 x zone, 实测列名 date
  UNION DISTINCT
  SELECT DISTINCT person_id, date AS dt
  FROM `{CDR}.activity_summary` WHERE steps >= 1  -- 日级; 弃用分钟级 steps_intraday 省扫描费
)
WHERE dt BETWEEN '{WIN_START}' AND '{WIN_END}'
""")
print(bq(f"SELECT COUNT(*) n, COUNT(DISTINCT person_id) p FROM `{OUT}.wear_days`"))

# %% Cell 5: 汇总核对 (贴回给 Claude)
print(bq(f"""
SELECT 'cohort_nights' t, COUNT(*) n, COUNT(DISTINCT person_id) p FROM `{OUT}.cohort_nights`
UNION ALL SELECT 'person_zip3', COUNT(*), COUNT(DISTINCT person_id) FROM `{OUT}.person_zip3`
UNION ALL SELECT 'person_covars', COUNT(*), COUNT(DISTINCT person_id) FROM `{OUT}.person_covars`
UNION ALL SELECT 'wear_days', COUNT(*), COUNT(DISTINCT person_id) FROM `{OUT}.wear_days`
"""))
