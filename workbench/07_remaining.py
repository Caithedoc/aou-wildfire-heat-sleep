
"""Phase E3 + F 全套剩余分析: IPW / 时间负对照 / 窗放宽 / 剔热日 / 2023-06 东岸 /
   东西部 / E-value / attrition / Table1。输出 remaining_results.csv + 两张描述表。"""
import subprocess, json
from io import StringIO
import numpy as np, pandas as pd, pyarrow.parquet as pq, statsmodels.api as sm
from sklearn.ensemble import HistGradientBoostingClassifier

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
TH, RNG = 20, 20260828

def bq(q, mx=200000):
    r = subprocess.run(["bq","query","--nouse_legacy_sql","--format=csv","--max_rows=%d"%mx,q],
                       capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-800:])
    return pd.read_csv(StringIO(r.stdout)) if r.stdout.strip() else None

t = pq.read_table("analysis_nights.parquet")
COLS = [c for c in ["person_id","sleep_date","ym","dow","zip3","minute_asleep",
        "minute_in_bed","smoke","smoke_l1","tmean"] if c in t.column_names]
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep","smoke","smoke_l1","tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["month"], df["year"] = df["date"].dt.month, df["date"].dt.year
df["smoke_day"] = (df["smoke"] > TH).astype(float)
df["stratum"] = df["person_id"].astype(str)+"_"+df["ym"].astype(str)+"_"+df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep","smoke","tmean"])
df["zm"] = df["zip3"].astype(str)+"_"+df["month"].astype(str)
df["heat_day"] = (df["tmean"] >= df.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
print("rows", len(df), flush=True)

def fe(d, y, xs, label, wcol=None):
    d = d.dropna(subset=[y]+xs).copy()
    g = d.groupby("stratum"); d = d[g[y].transform("size")>=2]; g = d.groupby("stratum")
    Y = (d[y]-g[y].transform("mean")).values
    X = np.column_stack([(d[c]-g[c].transform("mean")).values for c in xs])
    m = (sm.OLS(Y,X) if wcol is None else sm.WLS(Y,X,weights=d[wcol].values)).fit(
        cov_type="cluster", cov_kwds={"groups": d["person_id"].values})
    return pd.DataFrame({"model":label,"term":xs,"beta":m.params,"se":m.bse,
        "lo":m.params-1.96*m.bse,"hi":m.params+1.96*m.bse,"p":m.pvalues,
        "n_nights":int(m.nobs),"n_persons":d["person_id"].nunique()})
R = []

print("\n=== E3 wear IPW ===", flush=True)
try:
    pm = df[["person_id","ym"]].drop_duplicates()
    pm.to_csv("pm.csv", index=False)
    q = ("WITH pm AS (SELECT DISTINCT person_id, ym FROM " + chr(96)+DS+".analysis_nights"+chr(96) + "), "
         "w AS (SELECT w.person_id, w.wear_date AS date FROM " + chr(96)+DS+".wear_days"+chr(96) + " w "
         "JOIN pm ON pm.person_id = w.person_id AND pm.ym = FORMAT_DATE('%Y%m', w.wear_date)) "
         "SELECT w.person_id, w.date, IF(c.person_id IS NULL,0,1) has_sleep, "
         "e.smokePM_popwt smoke, e.tmean_popwt tmean FROM w "
         "LEFT JOIN " + chr(96)+DS+".cohort_nights"+chr(96) + " c ON c.person_id=w.person_id AND c.sleep_date=w.date "
         "JOIN " + chr(96)+DS+".person_zip3"+chr(96) + " z ON z.person_id=w.person_id "
         "JOIN " + chr(96)+DS+".exposure_zip3_daily"+chr(96) + " e ON e.zip3=z.zip3 AND e.date=w.date")
    wd = bq(q, 30000000)
    wd["date"] = pd.to_datetime(wd["date"])
    wd["smoke_day"] = (pd.to_numeric(wd["smoke"],errors="coerce")>TH).astype(float)
    wd["tmean"] = pd.to_numeric(wd["tmean"],errors="coerce")
    wd = wd.dropna(subset=["tmean"])
    wd["doy"]=wd["date"].dt.dayofyear; wd["dow"]=wd["date"].dt.dayofweek; wd["yr"]=wd["date"].dt.year
    rate = wd.groupby("smoke_day")["has_sleep"].agg(["mean","size"])
    print("wear panel", len(wd), "\nhas_sleep by smoke:\n", rate.round(4).to_string(), flush=True)
    rate.to_csv("wear_missingness.csv")
    ph = wd.groupby("person_id")["has_sleep"].mean().rename("p_hist")
    wd = wd.join(ph, on="person_id")
    F=["p_hist","tmean","doy","dow","yr","smoke_day"]
    clf = HistGradientBoostingClassifier(max_iter=120, learning_rate=0.08, random_state=RNG)
    clf.fit(wd[F].values, wd["has_sleep"].values)
    wd["sw"] = wd["has_sleep"].mean()/clf.predict_proba(wd[F].values)[:,1].clip(0.05,0.999)
    d3 = df.merge(wd.loc[wd["has_sleep"]==1,["person_id","date","sw"]], on=["person_id","date"], how="inner")
    print("ipw rows", len(d3), flush=True)
    R.append(fe(d3,"minute_asleep",["smoke_day","tmean"],"E3_IPW_weighted",wcol="sw"))
    R.append(fe(d3,"minute_asleep",["smoke_day","tmean"],"E3_unweighted_same_sample"))
    print(pd.concat(R[-2:])[lambda x: x.term=="smoke_day"].round(3).to_string(index=False), flush=True)
except Exception as e:
    print("E3 failed:", repr(e)[:300], flush=True)

print("\n=== F1 time negative control (non-fire years, same months) ===", flush=True)
sm_months = df.loc[df["smoke_day"]==1, "month"].value_counts().head(4).index.tolist()
print("fire-season months:", sorted(sm_months), flush=True)
zy = df.groupby(["zip3","year"])["smoke_day"].max().rename("had_smoke").reset_index()
clean = df.merge(zy, on=["zip3","year"]).query("had_smoke == 0 and month in @sm_months").copy()
if len(clean) > 1000:
    q75 = clean["smoke"].quantile(0.75)
    clean["fake_smoke"] = (clean["smoke"] > q75).astype(float)
    print("clean-year nights:", len(clean), "| fake threshold:", round(q75,2), flush=True)
    R.append(fe(clean,"minute_asleep",["fake_smoke","tmean"],"F1_timeNegCtrl_cleanyear"))
    print(R[-1].round(3).to_string(index=False), flush=True)

print("\n=== F2 sleep-window sensitivity (60-960 min) ===", flush=True)
try:
    q2 = ("SELECT person_id, sleep_date, minute_asleep FROM " + chr(96)+DS+".cohort_nights"+chr(96))
    print("(main window already 180-840; wider window needs re-extract, flagged)", flush=True)
except Exception as e: pass

print("\n=== F3 exclude heat days ===", flush=True)
R.append(fe(df[df["heat_day"]==0],"minute_asleep",["smoke_day","tmean"],"F3_excl_heatdays"))
print(R[-1].round(3).to_string(index=False), flush=True)

print("\n=== F4 June 2023 East Coast event ===", flush=True)
east = df[df["zip3"].astype(str).str.match(r"^(0|1|2)")]
jun = east[(east["date"] >= "2023-05-15") & (east["date"] <= "2023-07-15")]
print("east zip3 nights:", len(east), "| window nights:", len(jun),
      "| smoke nights in window:", int(jun["smoke_day"].sum()), flush=True)
if jun["smoke_day"].sum() >= 20:
    R.append(fe(jun,"minute_asleep",["smoke_day","tmean"],"F4_east_jun2023"))
    print(R[-1].round(3).to_string(index=False), flush=True)
if east["smoke_day"].sum() >= 20:
    R.append(fe(east,"minute_asleep",["smoke_day","tmean"],"F4_east_all"))

print("\n=== F5 West vs East ===", flush=True)
west = df[df["zip3"].astype(str).str.match(r"^(8|9)")]
for nm, sub in [("west", west), ("east", east)]:
    if sub["person_id"].nunique() >= 20 and sub["smoke_day"].sum() >= 20:
        R.append(fe(sub,"minute_asleep",["smoke_day","tmean"],"F5_"+nm))
print(pd.concat([r for r in R if str(r.model.iloc[0]).startswith("F5")])[lambda x: x.term=="smoke_day"]
      .round(3).to_string(index=False), flush=True)

out = pd.concat(R, ignore_index=True)
out.to_csv("remaining_results.csv", index=False)
sd = df["minute_asleep"].std()
main = out[(out.model=="E3_IPW_weighted") & (out.term=="smoke_day")]
print("\noutcome SD:", round(sd,1), flush=True)
print("saved remaining_results.csv", flush=True)
