"""13_ipw: compound IPW retry (2026-08-28). Full-panel bq CSV pull OOM/failed in 11;
retrain wear model on a 5% person-day sample + per-person p_hist from aggregated BQ,
then apply stabilized weights to analysis nights and refit compound interaction."""
import subprocess
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import statsmodels.api as sm
from sklearn.ensemble import HistGradientBoostingClassifier

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
THRESH = 20
RNG = 20260828

def bqcsv(sql, max_rows="30000000"):
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        f"--max_rows={max_rows}", sql], capture_output=True, text=True)
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[-300:]); return None
    return pd.read_csv(StringIO(r.stdout))

# per-person availability history (aggregated, 58K rows)
ph = bqcsv(f"""
WITH span AS (SELECT person_id, MIN(sleep_date) d0, MAX(sleep_date) d1
              FROM `{DS}.cohort_nights` GROUP BY person_id),
w AS (SELECT w.person_id, w.wear_date AS date
      FROM `{DS}.wear_days` w JOIN span s USING (person_id)
      WHERE w.wear_date BETWEEN s.d0 AND s.d1)
SELECT w.person_id,
       AVG(IF(c.person_id IS NULL, 0, 1)) AS p_hist,
       AVG(IF(c.person_id IS NULL, 0, 1)) AS obs_rate
FROM w LEFT JOIN `{DS}.cohort_nights` c
  ON c.person_id = w.person_id AND c.sleep_date = w.date
GROUP BY w.person_id""", max_rows="200000")
print("p_hist rows:", 0 if ph is None else len(ph))

# 5% sampled training panel
tr = bqcsv(f"""
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
WHERE MOD(ABS(FARM_FINGERPRINT(CONCAT(CAST(w.person_id AS STRING), CAST(w.date AS STRING)))), 20) = 0
""")
print("train rows:", 0 if tr is None else len(tr))

if ph is not None and tr is not None and len(tr):
    tr["date"] = pd.to_datetime(tr["date"])
    tr = tr.merge(ph[["person_id", "p_hist"]], on="person_id", how="left")
    for c in ("smoke", "tmean", "p_hist"):
        tr[c] = pd.to_numeric(tr[c], errors="coerce")
    tr = tr.dropna(subset=["tmean", "p_hist"])
    tr["smoke_day"] = (tr["smoke"] > THRESH).astype(float)
    tr["doy"] = tr["date"].dt.dayofyear
    tr["dw"] = tr["date"].dt.dayofweek
    tr["yr"] = tr["date"].dt.year
    F = ["p_hist", "tmean", "doy", "dw", "yr", "smoke_day"]
    clf = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.08,
                                         random_state=RNG)
    clf.fit(tr[F].values, tr["has_sleep"].values)
    marg = tr["has_sleep"].mean()
    print("sample has_sleep rate:", round(marg, 4))

    COLS = ["person_id", "sleep_date", "ym", "dow", "zip3",
            "minute_asleep", "smoke", "tmean"]
    t = pq.read_table("analysis_nights.parquet")
    df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
    for c in ["minute_asleep", "smoke", "tmean"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
    df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
    df["month"] = df["date"].dt.month
    df["smoke_day"] = (df["smoke"] > THRESH).astype(float)
    df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                     + "_" + df["dow"].astype(str))
    df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
    df["zm"] = df["zip3"].astype(str) + "_" + df["month"].astype(str)
    sz = df.groupby("stratum")["minute_asleep"].transform("size")
    df = df[sz >= 2].copy()
    df = df.merge(ph[["person_id", "p_hist"]], on="person_id", how="left")
    df["p_hist"] = pd.to_numeric(df["p_hist"], errors="coerce")
    df = df.dropna(subset=["p_hist"])
    df["doy"] = df["date"].dt.dayofyear
    df["dw"] = df["date"].dt.dayofweek
    df["yr"] = df["date"].dt.year
    df["p_obs"] = clf.predict_proba(df[F].values)[:, 1].clip(0.05, 0.999)
    df["sw"] = marg / df["p_obs"]
    thr = df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))
    df["hd"] = (df["tmean"] >= thr).astype(float)
    df["sxh"] = df["smoke_day"] * df["hd"]
    dd = df.copy()
    g = dd.groupby("stratum")
    dd = dd[g["minute_asleep"].transform("size") >= 2]
    g = dd.groupby("stratum")
    Y = (dd["minute_asleep"] - g["minute_asleep"].transform("mean")).values
    X = np.column_stack([(dd[c] - g[c].transform("mean")).values
                         for c in ["smoke_day", "hd", "sxh", "tmean"]])
    mw = sm.WLS(Y, X, weights=dd["sw"].values).fit(
        cov_type="cluster", cov_kwds={"groups": dd["person_id"].values})
    out = pd.DataFrame({"model": "cmp_ipw_sampled", "term": ["smoke_day", "hd", "sxh", "tmean"],
                        "beta": mw.params, "se": mw.bse,
                        "lo": mw.params - 1.96 * mw.bse, "hi": mw.params + 1.96 * mw.bse,
                        "p": mw.pvalues, "n": int(mw.nobs),
                        "n_both": int(dd["sxh"].sum())})
    print(out.round(4).to_string(index=False))
    out.to_csv("rev_cmp_ipw.csv", index=False)
    print("saved rev_cmp_ipw.csv")
