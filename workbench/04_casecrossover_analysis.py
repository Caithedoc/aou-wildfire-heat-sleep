"""Phase D: 挂接 + time-stratified case-crossover 主分析。
在 aou-wildfire-sleep 的 AoU Jupyter 里跑 (python3 04_casecrossover_analysis.py 或逐段贴)。

前置: wildfire_sleep 四表已建 (02); exposure csv 已上传 VM 并 bq load 为
      `wildfire_sleep.exposure_zip3_daily` (zip3 STRING, date DATE,
      smokePM_popwt FLOAT64, tmean_popwt FLOAT64)。

设计 (ANALYSIS_PLAN 预设):
- 层 = person x 年月 x 星期几 (time-stratified referent), 连续结局的
  case-crossover = 层内固定效应 OLS, 按人聚类稳健 SE。
- 主暴露 = 烟日 (zip3 日 smokePM > 20 µg/m3), 滞后 0-7 分布滞后;
  敏感性阈值 10/35 + 连续每 10 µg/m3。
- 温度: 层内已控季节, 再加当日 tmean 线性 + 平方 (敏感性: 不调 / 样条)。
- 效率: 只取"层内含 >=1 烟日"的 person-月 (其余层对烟系数无贡献)。
- 证伪: lead1 (明日烟) 不应预测今晚睡眠。
- 全程只导出系数表等聚合结果。
"""

import subprocess
from io import StringIO

import numpy as np
import pandas as pd

CDRP = "wb-clever-rhubarb-4688"
DS = f"{CDRP}.wildfire_sleep"
THRESH = 20  # 主定义; 敏感性 10 / 35


def bq(sql, quiet=False):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        "--max_rows=100000", sql], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-2000:])
    if quiet or not r.stdout.strip():
        return None
    return pd.read_csv(StringIO(r.stdout))


# %% Step 1: 暴露表加滞后列 (zip3 内按日期 LAG; lead1 做负对照)
bq(f"""
CREATE OR REPLACE TABLE `{DS}.exposure_lags` AS
SELECT zip3, date, smokePM_popwt, tmean_popwt,
  LAG(smokePM_popwt, 1) OVER w AS smoke_l1,
  LAG(smokePM_popwt, 2) OVER w AS smoke_l2,
  LAG(smokePM_popwt, 3) OVER w AS smoke_l3,
  LAG(smokePM_popwt, 4) OVER w AS smoke_l4,
  LAG(smokePM_popwt, 5) OVER w AS smoke_l5,
  LAG(smokePM_popwt, 6) OVER w AS smoke_l6,
  LAG(smokePM_popwt, 7) OVER w AS smoke_l7,
  LEAD(smokePM_popwt, 1) OVER w AS smoke_lead1,
  LAG(tmean_popwt, 1) OVER w AS tmean_l1
FROM `{DS}.exposure_zip3_daily`
WINDOW w AS (PARTITION BY zip3 ORDER BY date)
""", quiet=True)
print("exposure_lags built")

# %% Step 2: 挂接 + 取分析层 (person-月 内含烟日; 剔搬家不可行, zip 为快照)
bq(f"""
CREATE OR REPLACE TABLE `{DS}.analysis_nights` AS
WITH linked AS (
  SELECT c.person_id, c.sleep_date, c.minute_asleep, c.minute_awake,
         c.minute_deep, c.minute_rem, c.sleep_efficiency, c.src_id,
         z.zip3,
         e.smokePM_popwt AS smoke, e.smoke_l1, e.smoke_l2, e.smoke_l3,
         e.smoke_l4, e.smoke_l5, e.smoke_l6, e.smoke_l7, e.smoke_lead1,
         e.tmean_popwt AS tmean, e.tmean_l1,
         FORMAT_DATE('%Y%m', c.sleep_date) AS ym,
         EXTRACT(DAYOFWEEK FROM c.sleep_date) AS dow
  FROM `{DS}.cohort_nights` c
  JOIN `{DS}.person_zip3` z USING (person_id)
  JOIN `{DS}.exposure_lags` e
    ON e.zip3 = z.zip3 AND e.date = c.sleep_date
),
strata AS (
  SELECT person_id, ym
  FROM linked GROUP BY person_id, ym
  HAVING MAX(smoke) > {THRESH}
)
SELECT l.* FROM linked l JOIN strata USING (person_id, ym)
""", quiet=True)
print(bq(f"""SELECT COUNT(*) nights, COUNT(DISTINCT person_id) persons,
  COUNT(DISTINCT CONCAT(person_id,'-',ym)) person_months,
  COUNTIF(smoke > {THRESH}) smoke_nights FROM `{DS}.analysis_nights`"""))

# %% Step 3: 拉进 VM 内存 (数据不出 workspace), 层内去均值 + OLS 聚类 SE
import statsmodels.api as sm

import os

CACHE = "analysis_nights.parquet"   # VM 本地缓存, 重跑敏感性不必再下载 5 分钟
df = None
try:
    if os.path.exists(CACHE):
        df = pd.read_parquet(CACHE)
        print("loaded from cache")
    else:
        import pandas_gbq
        df = pandas_gbq.read_gbq(f"SELECT * FROM `{DS}.analysis_nights`",
                                 project_id=CDRP, progress_bar_type=None)
        df.to_parquet(CACHE, index=False)
except Exception:
    subprocess.run(["bash", "-c",
        f"bq extract --destination_format=CSV '{DS}.analysis_nights' "
        f"gs://$(wb resource resolve --id=ws_files 2>/dev/null)/an.csv"],
        capture_output=True, text=True)
if df is None:
    raise SystemExit("pandas_gbq 不可用, 改走 bq extract + 本地读 (见上一段)")
print("rows in memory:", len(df))

# pandas_gbq 返回可空扩展类型 (Int64/Float64/boolean), statsmodels 只吃 numpy;
# 数值列统一降成 float64, 缺失保持 NaN (2026-08-28 实测必需)。
NUMCOLS = ["minute_asleep", "minute_awake", "minute_deep", "minute_rem",
           "sleep_efficiency", "smoke", "smoke_l1", "smoke_l2", "smoke_l3",
           "smoke_l4", "smoke_l5", "smoke_l6", "smoke_l7", "smoke_lead1",
           "tmean", "tmean_l1"]
for c in NUMCOLS:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")

df["stratum"] = df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str)
for c in ["smoke", "smoke_l1", "smoke_l2", "smoke_l3", "smoke_l4",
          "smoke_l5", "smoke_l6", "smoke_l7", "smoke_lead1"]:
    df["b_" + c] = (df[c] > THRESH).astype(float)
df["tmean2"] = df["tmean"] ** 2


def fe_ols(dat, ycol, xcols, label):
    """层内去均值 OLS, 按 person 聚类 SE。返回系数表 DataFrame。"""
    d = dat.dropna(subset=[ycol] + xcols).copy()
    g = d.groupby("stratum")
    keep = g[ycol].transform("size") >= 2          # 层内至少 2 夜才有对照
    d = d[keep]
    g = d.groupby("stratum")
    Y = d[ycol] - g[ycol].transform("mean")
    X = pd.DataFrame({c: d[c] - g[c].transform("mean") for c in xcols})
    m = sm.OLS(Y.values, X.values).fit(cov_type="cluster",
                                       cov_kwds={"groups": d["person_id"].values})
    out = pd.DataFrame({"model": label, "outcome": ycol, "term": xcols,
                        "beta": m.params, "se": m.bse,
                        "lo": m.params - 1.96 * m.bse,
                        "hi": m.params + 1.96 * m.bse, "p": m.pvalues,
                        "n_nights": int(m.nobs),
                        "n_persons": d["person_id"].nunique()})
    print(out.round(3).to_string(index=False))
    return out

results = []
TADJ = ["tmean", "tmean2"]

# 主分析: lag0-1 累积 (b_smoke + b_smoke_l1), 结局 = 总睡眠时长
results.append(fe_ols(df, "minute_asleep", ["b_smoke", "b_smoke_l1"] + TADJ, "primary_lag01"))

# 分布滞后 0-7 (同模)
DL = ["b_smoke"] + [f"b_smoke_l{k}" for k in range(1, 8)]
results.append(fe_ols(df, "minute_asleep", DL + TADJ, "dlag07"))

# 次级结局
for y in ["minute_awake", "minute_deep", "minute_rem", "sleep_efficiency"]:
    results.append(fe_ols(df, y, ["b_smoke", "b_smoke_l1"] + TADJ, "secondary_lag01"))

# 连续暴露 (每 10 µg/m3)
df["smoke10"] = df["smoke"] / 10
df["smoke10_l1"] = df["smoke_l1"] / 10
results.append(fe_ols(df, "minute_asleep", ["smoke10", "smoke10_l1"] + TADJ, "continuous_per10"))

# 温度敏感性: 不调温度
results.append(fe_ols(df, "minute_asleep", ["b_smoke", "b_smoke_l1"], "no_temp_adj"))

# 负对照: 明日烟 (含真实滞后共调, lead 系数应近 0)
results.append(fe_ols(df, "minute_asleep",
                      ["b_smoke_lead1", "b_smoke", "b_smoke_l1"] + TADJ, "negctrl_lead1"))

allres = pd.concat(results, ignore_index=True)
allres.to_csv("phaseD_results.csv", index=False)
print("saved phaseD_results.csv (聚合系数表, 可导出)")

# %% Step 4: 阈值敏感性 10/35 与 3-14h 窗放宽 —— 改 THRESH / 清理窗后重跑本脚本。
