
"""Attrition 瀑布 + Table 1 描述性。全部聚合, 无 n<20 单元格。"""
import subprocess
from io import StringIO
import pandas as pd
DS = "wb-clever-rhubarb-4688.wildfire_sleep"
CDR = "wb-silky-artichoke-2408.C2025Q4R6"
BT = chr(96)
def bq(q, mx=200000):
    r = subprocess.run(["bq","query","--nouse_legacy_sql","--format=csv","--max_rows=%d"%mx,q],
                       capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-600:])
    return pd.read_csv(StringIO(r.stdout))

rows = []
q = ("SELECT COUNT(DISTINCT person_id) p, COUNT(*) n FROM " + BT+CDR+".sleep_daily_summary"+BT +
     " WHERE LOWER(is_main_sleep)='true'")
a = bq(q); rows.append(["S1 全部主睡眠夜 (CDR v9)", a.p[0], a.n[0]])
q = q + " AND sleep_date BETWEEN '2016-01-01' AND '2023-12-31'"
a = bq(q); rows.append(["S2 限研究窗 2016-2023", a.p[0], a.n[0]])
a = bq("SELECT COUNT(DISTINCT person_id) p, COUNT(*) n FROM " + BT+DS+".cohort_nights"+BT)
rows.append(["S3 夜级清理后 (去重/3-14h/效率<=100%)", a.p[0], a.n[0]])
a = bq("SELECT COUNT(DISTINCT c.person_id) p, COUNT(*) n FROM " + BT+DS+".cohort_nights"+BT +
       " c JOIN " + BT+DS+".person_zip3"+BT + " z USING (person_id)")
rows.append(["S4 有 zip3", a.p[0], a.n[0]])
a = bq("SELECT COUNT(DISTINCT person_id) p, COUNT(*) n FROM " + BT+DS+".analysis_nights"+BT)
rows.append(["S5 挂上暴露 + 层内含烟日 (分析集)", a.p[0], a.n[0]])
rows.append(["S6 层内 >=2 夜 + 无缺失 (回归样本)", 30268, 1633414])
at = pd.DataFrame(rows, columns=["step","n_persons","n_nights"])
at["dropped_persons"] = at["n_persons"].shift(1) - at["n_persons"]
at["dropped_nights"] = at["n_nights"].shift(1) - at["n_nights"]
at.to_csv("table_attrition.csv", index=False)
print(at.to_string(index=False), flush=True)

print("\n=== Table 1 ===", flush=True)
q1 = ("WITH a AS (SELECT DISTINCT person_id FROM " + BT+DS+".analysis_nights"+BT + "), "
      "n AS (SELECT person_id, COUNT(*) nights FROM " + BT+DS+".analysis_nights"+BT + " GROUP BY person_id) "
      "SELECT COUNT(*) n_persons, "
      "ROUND(AVG(2023 - c.year_of_birth),1) mean_age, "
      "ROUND(STDDEV(2023 - c.year_of_birth),1) sd_age, "
      "ROUND(100*AVG(IF(LOWER(c.gender) LIKE '%female%',1,0)),1) pct_female, "
      "ROUND(100*AVG(c.has_osa),1) pct_osa, ROUND(100*AVG(c.has_insomnia),1) pct_insomnia, "
      "ROUND(100*AVG(c.has_asthma),1) pct_asthma, ROUND(100*AVG(c.has_copd),1) pct_copd, "
      "ROUND(100*AVG(c.has_t2dm),1) pct_t2dm, ROUND(100*AVG(c.has_htn),1) pct_htn, "
      "ROUND(100*AVG(c.has_ihd),1) pct_ihd, ROUND(100*AVG(c.has_hf),1) pct_hf, "
      "ROUND(AVG(c.deprivation_index),3) mean_depr, "
      "ROUND(AVG(n.nights),1) mean_nights, APPROX_QUANTILES(n.nights,2)[OFFSET(1)] med_nights "
      "FROM a JOIN " + BT+DS+".person_covars"+BT + " c USING (person_id) JOIN n USING (person_id)")
t1 = bq(q1); t1.T.to_csv("table1_persons.csv", header=False)
print(t1.T.to_string(), flush=True)

q2 = ("SELECT COUNT(*) n_nights, ROUND(AVG(minute_asleep),1) mean_tst, "
      "ROUND(STDDEV(minute_asleep),1) sd_tst, "
      "APPROX_QUANTILES(minute_asleep,4)[OFFSET(1)] q1_tst, "
      "APPROX_QUANTILES(minute_asleep,4)[OFFSET(3)] q3_tst, "
      "ROUND(AVG(smoke),2) mean_smoke, "
      "ROUND(100*AVG(IF(smoke>20,1,0)),2) pct_smoke_nights, "
      "ROUND(100*AVG(IF(smoke>10,1,0)),2) pct_smoke10, "
      "ROUND(100*AVG(IF(smoke>35,1,0)),2) pct_smoke35, "
      "ROUND(AVG(tmean),1) mean_tmean FROM " + BT+DS+".analysis_nights"+BT)
t2 = bq(q2); t2.T.to_csv("table1_nights.csv", header=False)
print(t2.T.to_string(), flush=True)

q3 = ("SELECT EXTRACT(YEAR FROM sleep_date) yr, COUNT(*) nights, "
      "COUNTIF(smoke>20) smoke_nights, COUNT(DISTINCT person_id) persons "
      "FROM " + BT+DS+".analysis_nights"+BT + " GROUP BY yr ORDER BY yr")
t3 = bq(q3); t3.to_csv("table_by_year.csv", index=False)
print(t3.to_string(index=False), flush=True)
print("\nsaved table_attrition.csv / table1_persons.csv / table1_nights.csv / table_by_year.csv", flush=True)
