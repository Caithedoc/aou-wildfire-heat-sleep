"""12_reconcile: Table 3 定稿数据 + 18,168 来历对账 (2026-08-28)。
(a) E4 老口径复现: 阈值/计数都在 pre-filter 帧 (预期 both=18,168)
(b) 新口径 (回归帧) 四格模型完整输出 (带 CI) + 交互模型 -> rev_table3_final.csv
"""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import statsmodels.api as sm

THRESH = 20
COLS = ["person_id", "sleep_date", "ym", "dow", "zip3",
        "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
dfa = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    dfa[c] = pd.to_numeric(dfa[c], errors="coerce").astype("float64")
dfa["date"] = pd.to_datetime(dfa["sleep_date"].astype(str), errors="coerce")
dfa["month"] = dfa["date"].dt.month
dfa["smoke_day"] = (dfa["smoke"] > THRESH).astype(float)
dfa["stratum"] = (dfa["person_id"].astype(str) + "_" + dfa["ym"].astype(str)
                  + "_" + dfa["dow"].astype(str))
dfa = dfa.dropna(subset=["minute_asleep", "smoke", "tmean"])
dfa["zm"] = dfa["zip3"].astype(str) + "_" + dfa["month"].astype(str)
print("pre-filter frame:", len(dfa), "nights,", dfa["person_id"].nunique(), "persons")

def hf(d, q=0.90):
    thr = d.groupby("zm")["tmean"].transform(lambda s: s.quantile(q))
    return (d["tmean"] >= thr).astype(float)

# (a) E4 老口径: flags + counts on pre-filter frame
dfa["heat_e4"] = hf(dfa)
def cellcount(d, hcol):
    return dict(neither=int(((d[hcol]==0)&(d.smoke_day==0)).sum()),
                hot_only=int(((d[hcol]==1)&(d.smoke_day==0)).sum()),
                smoke_only=int(((d[hcol]==0)&(d.smoke_day==1)).sum()),
                both=int(((d[hcol]==1)&(d.smoke_day==1)).sum()), total=len(d))
print("E4-frame counts (pre-filter):", cellcount(dfa, "heat_e4"))
sz = dfa.groupby("stratum")["minute_asleep"].transform("size")
dfr = dfa[sz >= 2].copy()
print("E4 flags on model sample:", cellcount(dfr, "heat_e4"))

# (b) 新口径: flags recomputed on regression frame; full models with CI
dfr["heat90"] = hf(dfr)
print("consistent-frame counts:", cellcount(dfr, "heat90"))
def fe(d, y, xs, label):
    d = d.dropna(subset=[y]+xs).copy()
    g = d.groupby("stratum"); d = d[g[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y]-g[y].transform("mean")).values
    X = np.column_stack([(d[c]-g[c].transform("mean")).values for c in xs])
    m = sm.OLS(Y, X).fit(cov_type="cluster", cov_kwds={"groups": d["person_id"].values})
    out = pd.DataFrame({"model": label, "term": xs, "beta": m.params, "se": m.bse,
                        "lo": m.params-1.96*m.bse, "hi": m.params+1.96*m.bse,
                        "p": m.pvalues, "n": int(m.nobs)})
    print(out.round(4).to_string(index=False))
    return out

dfr["g_h"] = ((dfr.heat90==1)&(dfr.smoke_day==0)).astype(float)
dfr["g_s"] = ((dfr.heat90==0)&(dfr.smoke_day==1)).astype(float)
dfr["g_b"] = ((dfr.heat90==1)&(dfr.smoke_day==1)).astype(float)
dfr["sxh"] = dfr.smoke_day*dfr.heat90
o1 = fe(dfr, "minute_asleep", ["g_h","g_s","g_b"], "final_4cell")
o2 = fe(dfr, "minute_asleep", ["smoke_day","heat90","sxh","tmean"], "final_interaction")
pd.concat([o1,o2]).to_csv("rev_table3_final.csv", index=False)
print("saved rev_table3_final.csv")
