"""Phase A: CDR recon. 在 AoU Researcher Workbench 的 Jupyter (Python) 里逐 cell 跑。
回答 ANALYSIS_PLAN.md 的五个 recon 问题。只读,不落任何表。

用法: 复制进 Workbench notebook。所有表名/列名先由 discovery cell 实测,
不要信任何硬编码字段名 (包括本文件注释里提到的)。
"""

# %% Cell 1: 环境与 CDR 标识 (回答: tier + CDR 版本)
# 2026-06-30 起 AoU 已迁到 Researcher Workbench 2.0 (Verily Workbench)。
# 老平台的 WORKSPACE_CDR 环境变量在 2.0 可能不存在;2.0 里 CDR 是挂进
# workspace 的 referenced BigQuery dataset 资源,用 wb CLI 解析。
import os
import subprocess
import pandas as pd

CDR = os.environ.get("WORKSPACE_CDR")  # 老变量,能取到就直接用
if not CDR:
    # 列出 workspace 资源,找 All of Us Controlled Tier CDR 的资源名:
    listing = subprocess.run(["wb", "resource", "list"],
                             capture_output=True, text=True)
    print(listing.stdout or listing.stderr)
    RESOURCE_ID = "aou_cdr_ct_C2025Q4R6"  # 2026-08-28 已由 CLI 挂入 workspace
    # (指向 wb-silky-artichoke-2408.C2025Q4R6, Controlled Tier 2025Q4 切面, 133 表)
    if RESOURCE_ID:
        CDR = subprocess.run(["wb", "resolve", "--id", RESOURCE_ID],
                             capture_output=True, text=True).stdout.strip()
print("CDR =", CDR)  # 形如 "project_id.dataset_id"
# 判读: 资源/数据集名含 "controlled"/"ct" = Controlled Tier;
# 以 Data collection 页面标注 (All of Us Controlled Tier, CDR v8) 为最终答案。


def gbq(sql: str) -> pd.DataFrame:
    return pd.read_gbq(sql, dialect="standard", progress_bar_type=None)


# %% Cell 2: Fitbit 相关表都有哪些 (discovery,不猜表名)
tables = gbq(f"""
SELECT table_name
FROM `{CDR}.INFORMATION_SCHEMA.TABLES`
WHERE LOWER(table_name) LIKE '%sleep%'
   OR LOWER(table_name) LIKE '%fitbit%'
   OR LOWER(table_name) LIKE '%activity%'
   OR LOWER(table_name) LIKE '%heart%'
   OR LOWER(table_name) LIKE '%steps%'
   OR LOWER(table_name) LIKE '%zip%'
ORDER BY table_name
""")
print(tables.to_string())

# %% Cell 3: sleep_daily_summary 的列 (discovery)
SLEEP_TABLE = "sleep_daily_summary"  # 若 Cell 2 显示不同名字,在此改
cols = gbq(f"""
SELECT column_name, data_type
FROM `{CDR}.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name = '{SLEEP_TABLE}'
ORDER BY ordinal_position
""")
print(cols.to_string())

# %% Cell 4: 日期范围 + 人数 + 人夜数 (回答 recon 问题 2 与 3)
# 列名按 Cell 3 实测调整;常见为 person_id / sleep_date / is_main_sleep / minute_asleep
summary = gbq(f"""
SELECT
  MIN(sleep_date) AS min_date,
  MAX(sleep_date) AS max_date,
  COUNT(DISTINCT person_id) AS n_persons,
  COUNT(*) AS n_person_nights
FROM `{CDR}.{SLEEP_TABLE}`
""")
print(summary.to_string())
print("\n关键判读: max_date 是否 >= 2023-06-30 (东岸烟事件窗)?")

# %% Cell 5: 逐年人夜数 (看 BYOD 段 vs 设备发放段的量变)
by_year = gbq(f"""
SELECT EXTRACT(YEAR FROM sleep_date) AS yr,
       COUNT(DISTINCT person_id) AS n_persons,
       COUNT(*) AS n_nights
FROM `{CDR}.{SLEEP_TABLE}`
GROUP BY yr ORDER BY yr
""")
print(by_year.to_string())

# %% Cell 6: zip3 可得性 (回答 recon 问题 4)
# Controlled Tier 的 zip3 常见于 observation 表 (PII_STATE / zip 相关 concept) 或专表
# (如 zip3_ses_map 关联)。先 discovery:
zip_tables = gbq(f"""
SELECT table_name FROM `{CDR}.INFORMATION_SCHEMA.TABLES`
WHERE LOWER(table_name) LIKE '%zip%' OR LOWER(table_name) LIKE '%ses%'
""")
print(zip_tables.to_string())
# 找到 zip3 来源表后,统计: 有 Fitbit 睡眠且有 zip3 的人数,按 zip3 前缀分布
# (西部烟区 8xx/9xx,东北 0xx/1xx),此处按实测表名补 SQL。

# %% Cell 7: 佩戴天数分布 (回答 recon 问题 5,给 IPW 打底)
wear = gbq(f"""
SELECT person_id, COUNT(DISTINCT sleep_date) AS n_nights
FROM `{CDR}.{SLEEP_TABLE}`
GROUP BY person_id
""")
print(wear["n_nights"].describe())
print("\n>=180 夜人数:", int((wear["n_nights"] >= 180).sum()))
print(">=365 夜人数:", int((wear["n_nights"] >= 365).sum()))

# %% Cell 8: recon 结论汇总 (手填后贴回给 Claude 排 Phase B)
RECON = {
    "tier": "",            # Registered / Controlled
    "cdr_version": "",
    "sleep_date_range": "",
    "covers_2023_june_east_coast": None,  # True/False, 决定头号自然实验
    "n_persons_sleep": None,
    "n_person_nights": None,
    "zip3_source_table": "",
    "n_persons_ge180_nights": None,
}
print(RECON)
