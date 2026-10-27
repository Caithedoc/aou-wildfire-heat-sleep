"""Phase F: 阈值敏感性 10/35 + 设备来源分层。绕开 dbdate: arrow 列直接转 numpy。"""
import pandas as pd, numpy as np, statsmodels.api as sm, pyarrow.parquet as pq
COLS = ["person_id","ym","dow","src_id","minute_asleep","smoke","smoke_l1","tmean"]
t = pq.read_table("analysis_nights.parquet", columns=COLS)
df0 = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep","smoke","smoke_l1","tmean"]:
    df0[c] = pd.to_numeric(df0[c], errors="coerce").astype("float64")
df0["stratum"] = df0["person_id"].astype(str)+"_"+df0["ym"].astype(str)+"_"+df0["dow"].astype(str)
df0["tmean2"] = df0["tmean"]**2
print("cache rows", len(df0))

def fe(dat, y, xs, label):
    d = dat.dropna(subset=[y]+xs).copy()
    g = d.groupby("stratum")
    d = d[g[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    m = sm.OLS(Y, X).fit(cov_type="cluster", cov_kwds={"groups": d["person_id"].values})
    return pd.DataFrame({"model":label,"term":xs,"beta":m.params,"se":m.bse,
                         "lo":m.params-1.96*m.bse,"hi":m.params+1.96*m.bse,"p":m.pvalues,
                         "n_nights":int(m.nobs),"n_persons":d["person_id"].nunique()})

res = []
for TH in [10, 35]:
    d = df0[df0.groupby(["person_id","ym"])["smoke"].transform("max") > TH].copy()
    d["b0"] = (d["smoke"] > TH).astype(float); d["b1"] = (d["smoke_l1"] > TH).astype(float)
    res.append(fe(d, "minute_asleep", ["b0","b1","tmean","tmean2"], "thresh_%d" % TH))

d = df0.copy()
d["b0"] = (d["smoke"] > 20).astype(float); d["b1"] = (d["smoke_l1"] > 20).astype(float)
for src, g in d.groupby("src_id"):
    if g["person_id"].nunique() < 500: continue
    res.append(fe(g, "minute_asleep", ["b0","b1","tmean","tmean2"], "src_" + str(src)[:20]))

out = pd.concat(res, ignore_index=True)
out.to_csv("phaseF_results.csv", index=False)
print(out[out["term"].isin(["b0","b1"])].round(3).to_string(index=False))