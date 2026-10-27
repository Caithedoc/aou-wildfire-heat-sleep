
"""收尾: F1b 时间负对照(改定义) / F2 睡眠窗放宽 / E-value / 2023 剔除敏感性。"""
import subprocess
from io import StringIO
import numpy as np, pandas as pd, pyarrow.parquet as pq, statsmodels.api as sm
DS="wb-clever-rhubarb-4688.wildfire_sleep"; CDR="wb-silky-artichoke-2408.C2025Q4R6"
BT=chr(96); TH=20
def bq(q, mx=5000000):
    r=subprocess.run(["bq","query","--nouse_legacy_sql","--format=csv","--max_rows=%d"%mx,q],
                     capture_output=True,text=True)
    if r.returncode: raise RuntimeError(r.stderr[-500:])
    return pd.read_csv(StringIO(r.stdout)) if r.stdout.strip() else None

t=pq.read_table("analysis_nights.parquet")
C=[c for c in ["person_id","sleep_date","ym","dow","zip3","minute_asleep","smoke","tmean"] if c in t.column_names]
df=pd.DataFrame({c:t.column(c).to_numpy(zero_copy_only=False) for c in C})
for c in ["minute_asleep","smoke","tmean"]: df[c]=pd.to_numeric(df[c],errors="coerce").astype("float64")
df["date"]=pd.to_datetime(df["sleep_date"].astype(str),errors="coerce")
df["month"],df["year"]=df["date"].dt.month,df["date"].dt.year
df["smoke_day"]=(df["smoke"]>TH).astype(float)
df["stratum"]=df["person_id"].astype(str)+"_"+df["ym"].astype(str)+"_"+df["dow"].astype(str)
df=df.dropna(subset=["minute_asleep","smoke","tmean"])

def fe(d,y,xs,label):
    d=d.dropna(subset=[y]+xs).copy(); g=d.groupby("stratum")
    d=d[g[y].transform("size")>=2]; g=d.groupby("stratum")
    Y=(d[y]-g[y].transform("mean")).values
    X=np.column_stack([(d[c]-g[c].transform("mean")).values for c in xs])
    m=sm.OLS(Y,X).fit(cov_type="cluster",cov_kwds={"groups":d["person_id"].values})
    return pd.DataFrame({"model":label,"term":xs,"beta":m.params,"se":m.bse,
        "lo":m.params-1.96*m.bse,"hi":m.params+1.96*m.bse,"p":m.pvalues,
        "n_nights":int(m.nobs),"n_persons":d["person_id"].nunique()})
R=[]

print("=== F1b time negative control (low-smoke zip3-years, fire months) ===",flush=True)
mx=df.groupby(["zip3","year"])["smoke"].max().rename("mx").reset_index()
low=df.merge(mx,on=["zip3","year"]).query("mx < 10 and month in [6,7,8,9]").copy()
print("low-smoke-year fire-month nights:",len(low),"persons:",low["person_id"].nunique(),flush=True)
if len(low)>5000:
    thr=low["smoke"].quantile(0.90)
    low["fake"]=(low["smoke"]>thr).astype(float)
    print("fake threshold (p90 of low-smoke):",round(thr,2),flush=True)
    R.append(fe(low,"minute_asleep",["fake","tmean"],"F1b_timeNegCtrl"))
    print(R[-1].round(3).to_string(index=False),flush=True)

print("\n=== F2 wider sleep window 60-960 min ===",flush=True)
q=("WITH ms AS (SELECT person_id, sleep_date, minute_asleep, minute_in_bed, "
   "ROW_NUMBER() OVER (PARTITION BY person_id, sleep_date ORDER BY minute_asleep DESC) rk "
   "FROM "+BT+CDR+".sleep_daily_summary"+BT+" WHERE LOWER(is_main_sleep)='true' "
   "AND sleep_date BETWEEN '2016-01-01' AND '2023-12-31'), "
   "cn AS (SELECT person_id, sleep_date, minute_asleep FROM ms WHERE rk=1 "
   "AND minute_asleep BETWEEN 60 AND 960 AND minute_in_bed>=minute_asleep AND minute_in_bed>0), "
   "j AS (SELECT c.person_id, c.sleep_date, c.minute_asleep, e.smokePM_popwt smoke, e.tmean_popwt tmean, "
   "FORMAT_DATE('%Y%m', c.sleep_date) ym, EXTRACT(DAYOFWEEK FROM c.sleep_date) dow "
   "FROM cn c JOIN "+BT+DS+".person_zip3"+BT+" z USING (person_id) "
   "JOIN "+BT+DS+".exposure_zip3_daily"+BT+" e ON e.zip3=z.zip3 AND e.date=c.sleep_date), "
   "st AS (SELECT person_id, ym FROM j GROUP BY person_id, ym HAVING MAX(smoke)>20) "
   "SELECT j.* FROM j JOIN st USING (person_id, ym)")
w=bq(q, 5000000)
if w is not None and len(w):
    for c in ["minute_asleep","smoke","tmean"]: w[c]=pd.to_numeric(w[c],errors="coerce")
    w["smoke_day"]=(w["smoke"]>TH).astype(float)
    w["stratum"]=w["person_id"].astype(str)+"_"+w["ym"].astype(str)+"_"+w["dow"].astype(str)
    w=w.dropna(subset=["minute_asleep","smoke","tmean"])
    print("wide-window nights:",len(w),flush=True)
    R.append(fe(w,"minute_asleep",["smoke_day","tmean"],"F2_window_60_960"))
    print(R[-1].round(3).to_string(index=False),flush=True)

print("\n=== F6 exclude 2023 (event-dominance check) ===",flush=True)
R.append(fe(df[df["year"]!=2023],"minute_asleep",["smoke_day","tmean"],"F6_excl_2023"))
print(R[-1].round(3).to_string(index=False),flush=True)

print("\n=== E-values ===",flush=True)
SD=df["minute_asleep"].std()
def evalue(beta, sd=SD):
    d=abs(beta)/sd; rr=np.exp(0.91*d)
    return rr+np.sqrt(rr*(rr-1))
ev=pd.DataFrame({"estimate":["main_smoke_-0.58","compound_both_-4.02","cardioresp_-1.48","heat_per1C_-0.45"],
                 "beta":[-0.58,-4.02,-1.48,-0.45]})
ev["E_value"]=[round(evalue(b),3) for b in ev["beta"]]
ev["outcome_SD"]=round(SD,1)
print(ev.to_string(index=False),flush=True)
ev.to_csv("evalues.csv",index=False)

pd.concat(R,ignore_index=True).to_csv("final_sens.csv",index=False)
print("\nsaved final_sens.csv + evalues.csv",flush=True)
