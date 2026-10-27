
"""F1c 时间负对照 (正确版): 从完整 cohort_nights 取低烟 zip3-年的火季月, 分析集本身
   只含有烟月份, 从中找无烟年必然为空 —— 这是前两版返回 0 夜的原因。"""
import subprocess
from io import StringIO
import numpy as np, pandas as pd, statsmodels.api as sm
DS="wb-clever-rhubarb-4688.wildfire_sleep"; BT=chr(96)
def bq(q,mx=8000000):
    r=subprocess.run(["bq","query","--nouse_legacy_sql","--format=csv","--max_rows=%d"%mx,q],
                     capture_output=True,text=True)
    if r.returncode: raise RuntimeError(r.stderr[-500:])
    return pd.read_csv(StringIO(r.stdout)) if r.stdout.strip() else None

q=("WITH zy AS (SELECT zip3, EXTRACT(YEAR FROM date) yr, MAX(smokePM_popwt) mx "
   "FROM "+BT+DS+".exposure_zip3_daily"+BT+" GROUP BY zip3, yr HAVING mx < 5), "
   "j AS (SELECT c.person_id, c.sleep_date, c.minute_asleep, e.smokePM_popwt smoke, "
   "e.tmean_popwt tmean, FORMAT_DATE('%Y%m', c.sleep_date) ym, "
   "EXTRACT(DAYOFWEEK FROM c.sleep_date) dow, EXTRACT(MONTH FROM c.sleep_date) mo "
   "FROM "+BT+DS+".cohort_nights"+BT+" c "
   "JOIN "+BT+DS+".person_zip3"+BT+" z USING (person_id) "
   "JOIN "+BT+DS+".exposure_zip3_daily"+BT+" e ON e.zip3=z.zip3 AND e.date=c.sleep_date "
   "JOIN zy ON zy.zip3=z.zip3 AND zy.yr=EXTRACT(YEAR FROM c.sleep_date)) "
   "SELECT * FROM j WHERE mo BETWEEN 6 AND 9")
d=bq(q)
print("low-smoke-year fire-month nights:", 0 if d is None else len(d), flush=True)
if d is not None and len(d)>5000:
    for c in ["minute_asleep","smoke","tmean"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    d["stratum"]=d["person_id"].astype(str)+"_"+d["ym"].astype(str)+"_"+d["dow"].astype(str)
    d=d.dropna(subset=["minute_asleep","smoke","tmean"])
    thr=d["smoke"].quantile(0.95)
    d["fake"]=(d["smoke"]>thr).astype(float)
    print("persons:",d["person_id"].nunique()," fake-exposure threshold (p95):",round(thr,3),
          " fake-exposed nights:",int(d["fake"].sum()),flush=True)
    g=d.groupby("stratum"); d=d[g["minute_asleep"].transform("size")>=2]; g=d.groupby("stratum")
    Y=(d["minute_asleep"]-g["minute_asleep"].transform("mean")).values
    X=np.column_stack([(d[c]-g[c].transform("mean")).values for c in ["fake","tmean"]])
    m=sm.OLS(Y,X).fit(cov_type="cluster",cov_kwds={"groups":d["person_id"].values})
    res=pd.DataFrame({"model":"F1c_timeNegCtrl_lowSmokeYears","term":["fake","tmean"],
        "beta":m.params,"se":m.bse,"lo":m.params-1.96*m.bse,"hi":m.params+1.96*m.bse,
        "p":m.pvalues,"n_nights":int(m.nobs),"n_persons":d["person_id"].nunique()})
    print(res.round(3).to_string(index=False),flush=True)
    res.to_csv("timenegctrl.csv",index=False)
print("done",flush=True)
