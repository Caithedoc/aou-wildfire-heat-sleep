
"""Phase E: E1 baseline / E2 HTE+GATES / E3 wear-IPW / E4 smoke x heat joint."""
import subprocess
from io import StringIO
import numpy as np, pandas as pd, pyarrow.parquet as pq, statsmodels.api as sm
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.model_selection import KFold

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
TH, RNG = 20, 20260828

t = pq.read_table("analysis_nights.parquet")
COLS = [c for c in ["person_id","sleep_date","ym","dow","src_id","zip3",
        "minute_asleep","minute_awake","sleep_efficiency","smoke","smoke_l1","tmean"]
        if c in t.column_names]
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep","minute_awake","sleep_efficiency","smoke","smoke_l1","tmean"]:
    if c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["month"], df["doy"], df["year"] = df["date"].dt.month, df["date"].dt.dayofyear, df["date"].dt.year
df["smoke_day"] = (df["smoke"] > TH).astype(float)
df["stratum"] = df["person_id"].astype(str)+"_"+df["ym"].astype(str)+"_"+df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep","smoke","tmean"])
print("rows:", len(df), "persons:", df["person_id"].nunique(), flush=True)

def fe(d, y, xs, label, wcol=None):
    d = d.dropna(subset=[y]+xs).copy()
    g = d.groupby("stratum"); d = d[g[y].transform("size") >= 2]; g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    m = (sm.OLS(Y,X) if wcol is None else sm.WLS(Y,X,weights=d[wcol].values)).fit(
        cov_type="cluster", cov_kwds={"groups": d["person_id"].values})
    return pd.DataFrame({"model":label,"term":xs,"beta":m.params,"se":m.bse,
        "lo":m.params-1.96*m.bse,"hi":m.params+1.96*m.bse,"p":m.pvalues,
        "n_nights":int(m.nobs),"n_persons":d["person_id"].nunique()})

R = []

print("\n=== E1 counterfactual baseline ===", flush=True)
hist = (df[df["smoke"]<=5].groupby("person_id")["minute_asleep"].agg(["mean","std","count"])
        .rename(columns={"mean":"hist_mean","std":"hist_sd","count":"hist_n"}))
d1 = df.join(hist, on="person_id").dropna(subset=["hist_mean"])
d1 = d1[d1["hist_n"] >= 20]
F = ["hist_mean","hist_sd","tmean","doy","dow","year"]
tr = d1[d1["smoke"]<=5]
mdl = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, random_state=RNG)
mdl.fit(tr[F].values, tr["minute_asleep"].values)
d1["resid"] = d1["minute_asleep"] - mdl.predict(d1[F].values)
r1 = fe(d1, "resid", ["smoke_day","tmean"], "E1_residual_vs_expected")
print(r1.round(3).to_string(index=False), flush=True); R.append(r1)

print("\n=== E2 HTE + GATES ===", flush=True)
q = ("SELECT person_id, year_of_birth, deprivation_index, median_income, has_osa, "
     "has_insomnia, has_asthma, has_copd, has_t2dm, has_htn, has_ihd, has_hf FROM " + chr(96)+DS+".person_covars"+chr(96))
rr = subprocess.run(["bq","query","--nouse_legacy_sql","--format=csv","--max_rows=200000",q],
                    capture_output=True, text=True)
cov = pd.read_csv(StringIO(rr.stdout)) if rr.returncode == 0 else None
if cov is not None:
    d2 = df.merge(cov, on="person_id", how="inner").join(hist, on="person_id")
    d2["age"] = d2["year"] - d2["year_of_birth"]
    MODS = ["age","deprivation_index","median_income","hist_mean","hist_sd","has_osa",
            "has_insomnia","has_asthma","has_copd","has_t2dm","has_htn","has_ihd","has_hf"]
    d2 = d2.dropna(subset=["age","hist_mean"]).copy()
    for m in MODS: d2[m] = pd.to_numeric(d2[m], errors="coerce")
    d2[MODS] = d2[MODS].fillna(d2[MODS].median())
    g = d2.groupby("stratum"); d2 = d2[g["minute_asleep"].transform("size")>=2]; g = d2.groupby("stratum")
    d2["y_t"] = d2["minute_asleep"] - g["minute_asleep"].transform("mean")
    d2["w_t"] = d2["smoke_day"] - g["smoke_day"].transform("mean")
    f2 = d2[d2["w_t"].abs() > 1e-8].copy()
    f2["pseudo"] = f2["y_t"]/f2["w_t"]; f2["wt"] = f2["w_t"]**2
    tau = np.zeros(len(f2)); uniq = pd.unique(f2["person_id"].values)
    for tr_i, te_i in KFold(3, shuffle=True, random_state=RNG).split(uniq):
        mtr = f2["person_id"].isin(set(uniq[tr_i])).values
        mte = f2["person_id"].isin(set(uniq[te_i])).values
        ff = HistGradientBoostingRegressor(max_iter=150, learning_rate=0.08, random_state=RNG)
        ff.fit(f2.loc[mtr,MODS].values, f2.loc[mtr,"pseudo"].values, sample_weight=f2.loc[mtr,"wt"].values)
        tau[mte] = ff.predict(f2.loc[mte,MODS].values)
    f2["tau_hat"] = tau
    f2["grp"] = pd.qcut(f2["tau_hat"], 5, labels=False, duplicates="drop")
    gts = []
    for qq, gq in f2.groupby("grp"):
        sub = d2[d2["person_id"].isin(set(gq["person_id"]))]
        if sub["person_id"].nunique() < 20: continue
        x = fe(sub, "minute_asleep", ["smoke_day","tmean"], "E2_GATES_q%d" % (int(qq)+1))
        x["pred_tau"] = gq["tau_hat"].mean(); gts.append(x)
    if gts:
        gt = pd.concat(gts, ignore_index=True)
        print(gt[gt["term"]=="smoke_day"].round(3).to_string(index=False), flush=True); R.append(gt)
    SUBS = {"age_ge65": d2["age"]>=65, "age_lt65": d2["age"]<65,
            "osa": d2["has_osa"]==1, "no_osa": d2["has_osa"]==0,
            "cardioresp": d2[["has_asthma","has_copd","has_ihd","has_hf"]].sum(axis=1)>0,
            "deprived_top25": d2["deprivation_index"]>=d2["deprivation_index"].quantile(0.75)}
    sg = []
    for nm, mk in SUBS.items():
        sub = d2[mk]
        if sub["person_id"].nunique() < 20: continue
        sg.append(fe(sub, "minute_asleep", ["smoke_day","tmean"], "E2_sub_"+nm))
    if sg:
        s = pd.concat(sg, ignore_index=True)
        print(s[s["term"]=="smoke_day"].round(3).to_string(index=False), flush=True); R.append(s)

print("\n=== E4 smoke x heat joint ===", flush=True)
d4 = df.copy()
d4["zm"] = d4["zip3"].astype(str)+"_"+d4["month"].astype(str)
d4["heat_day"] = (d4["tmean"] >= d4.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
d4["sxh"] = d4["smoke_day"]*d4["heat_day"]
print("heat share:", round(d4["heat_day"].mean(),3), "both share:", round(d4["sxh"].mean(),4), flush=True)
x = fe(d4, "minute_asleep", ["smoke_day","heat_day","sxh","tmean"], "E4_interaction")
print(x.round(3).to_string(index=False), flush=True); R.append(x)
d4["g_heat"] = ((d4["heat_day"]==1)&(d4["smoke_day"]==0)).astype(float)
d4["g_smoke"] = ((d4["heat_day"]==0)&(d4["smoke_day"]==1)).astype(float)
d4["g_both"] = ((d4["heat_day"]==1)&(d4["smoke_day"]==1)).astype(float)
print("both-cell nights:", int(d4["g_both"].sum()), flush=True)
if d4["g_both"].sum() >= 20:
    x = fe(d4, "minute_asleep", ["g_heat","g_smoke","g_both"], "E4_joint_4cell")
    print(x.round(3).to_string(index=False), flush=True); R.append(x)

pd.concat(R, ignore_index=True).to_csv("phaseE_partial.csv", index=False)
print("\nsaved phaseE_partial.csv", flush=True)
