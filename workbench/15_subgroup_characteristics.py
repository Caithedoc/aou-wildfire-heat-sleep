"""15_subgroup_characteristics: 亚组特征表 + eTable 6 缺失的两行边际估计 (2026-09-24).

在 AoU Jupyter VM 里跑: `nohup python -u 15_subgroup_characteristics.py > rev5.log 2>&1 &`
需要同目录: analysis_nights.parquet (缓存, 同 11/12/13); person_covars 经 bq 拉取 (同 11), 失败时读 person_covars_local.csv。

亚组定义与 06_phaseE_ml.py / 11_revision.py 完全一致 (否则与 eTable 6 对不上):
  cardioresp     人级: has_asthma/has_copd/has_ihd/has_hf 任一 = 1        (事后定义)
  osa            人级: has_osa = 1
  age65          夜级: (year - year_of_birth) >= 65   → 同一人可能跨组 (Table 4 两组人数之和 30,640 > 30,268 即由此而来)
  deprived_top25 夜级: deprivation_index >= 分析夜分布的 P75 (阈值在夜上取, 同 11_revision.py:260)

输出 (聚合; 1-20 计数置 NaN):
  rev5_subgroup_characteristics.csv   每个亚组水平: 人数, 夜数, 首夜年龄 mean/SD, 女性 n(%), 地区夜数分布 %,
                                      烟雾夜 / 高温夜 / 复合夜 数量与占比, 烟雾夜平均 PM2.5
  rev5_subgroup_marginal.csv          8 个水平的边际烟雾模型 (smoke_day + tmean, person 聚类, 同 06 E2_subgroup):
                                      用于补 eTable 6 Panel A 的 "No cardiorespiratory disease" 与 "Less deprived quartiles" 两行,
                                      其余 6 行应与 eTable 6 现值一致 (作为核对)
  rev5_meta.csv
"""

import os, subprocess, time
from io import StringIO
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
THRESH = 20
MIN_CELL = 20
T0 = time.time()


def log(*a):
    print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)


def bqcsv(sql, max_rows="100000000"):
    try:
        r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                            f"--max_rows={max_rows}", sql], capture_output=True, text=True)
    except Exception as e:
        print("BQ UNAVAILABLE:", e); return None
    if r.returncode != 0:
        print("BQ FAIL:", r.stderr[:500]); return None
    return pd.read_csv(StringIO(r.stdout))


def sup(x):
    """1-20 计数抑制。"""
    return np.nan if (0 < x < MIN_CELL) else x


# ---------------------------------------------------------------- 读缓存 (同 13)
COLS = ["person_id", "sleep_date", "ym", "dow", "zip3", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["year"] = df["date"].dt.year
df["month"] = df["date"].dt.month
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
df["stratum"] = df["person_id"].astype(str) + "_" + df["ym"].astype(str) + "_" + df["dow"].astype(str)
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
df["zm"] = df["zip3"] + "_" + df["month"].astype(str)
df = df[df.groupby("stratum")["minute_asleep"].transform("size") >= 2].copy()
log("rows:", len(df), "| persons:", df["person_id"].nunique())

cov = bqcsv(f"SELECT person_id, year_of_birth, gender, deprivation_index, has_osa, "
            f"has_asthma, has_copd, has_ihd, has_hf FROM `{DS}.person_covars`", max_rows="200000")
if cov is None and os.path.exists("person_covars_local.csv"):
    cov = pd.read_csv("person_covars_local.csv")
assert cov is not None, "person_covars 不可用"
d = df.merge(cov, on="person_id", how="inner")
log("merged rows:", len(d), "| persons:", d["person_id"].nunique(), "(应与上行一致)")

# ---------------------------------------------------------------- 暴露与亚组 (同 06/11/13)
d["sd"] = (d["smoke"] > THRESH).astype(float)
d["hd"] = (d["tmean"] >= d.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))).astype(float)
d["joint"] = d["sd"] * d["hd"]
d["age"] = d["year"] - d["year_of_birth"]
d["cardioresp"] = (d[["has_asthma", "has_copd", "has_ihd", "has_hf"]].sum(axis=1) > 0)
d["osa"] = d["has_osa"] == 1
d["age65"] = d["age"] >= 65
thr75 = d["deprivation_index"].quantile(0.75)
d["deprived_top25"] = d["deprivation_index"] >= thr75
f = d["zip3"].str[0].astype(int)
d["region"] = np.select([f <= 2, f <= 7], ["East", "Central"], "West")
d["female"] = d["gender"] == "Female"
log(f"ADI P75 阈值 (夜分布) = {thr75:.3f}")

LEVELS = [  # (维度, 水平标签, 掩码)
    ("Cardiorespiratory disease", "Yes", d["cardioresp"]),
    ("Cardiorespiratory disease", "No", ~d["cardioresp"]),
    ("Obstructive sleep apnea", "Yes", d["osa"]),
    ("Obstructive sleep apnea", "No", ~d["osa"]),
    ("Age", "<65 years", ~d["age65"]),
    ("Age", "65 years or older", d["age65"]),
    ("Neighborhood deprivation", "Most deprived quartile", d["deprived_top25"]),
    ("Neighborhood deprivation", "Less deprived quartiles", ~d["deprived_top25"]),
]

# ---------------------------------------------------------------- 特征表
rows = []
for dim, lab, m in LEVELS + [("All", "All", pd.Series(True, index=d.index))]:
    s = d[m]
    first = s.sort_values("date").groupby("person_id").first()
    n_p = s["person_id"].nunique(); n_n = len(s)
    reg = s["region"].value_counts(normalize=True) * 100
    n_sd, n_hd, n_j = int(s["sd"].sum()), int(s["hd"].sum()), int(s["joint"].sum())
    rows.append({
        "dimension": dim, "level": lab,
        "participants": sup(n_p), "nights": sup(n_n),
        "age_first_night_mean": first["age"].mean(), "age_first_night_sd": first["age"].std(),
        "female_n": sup(int(first["female"].sum())), "female_pct": first["female"].mean() * 100,
        "region_east_pct": reg.get("East", 0), "region_central_pct": reg.get("Central", 0),
        "region_west_pct": reg.get("West", 0),
        "smoke_nights": sup(n_sd), "smoke_nights_pct": n_sd / n_n * 100,
        "hot_nights": sup(n_hd), "hot_nights_pct": n_hd / n_n * 100,
        "joint_nights": sup(n_j), "joint_nights_pct": n_j / n_n * 100,
        "smoke_pm_on_smoke_nights_mean": s.loc[s["sd"] == 1, "smoke"].mean(),
        "nights_per_participant_median": s.groupby("person_id").size().median(),
    })
ch = pd.DataFrame(rows)
ch.to_csv("rev5_subgroup_characteristics.csv", index=False)
log("=== 亚组特征 ===\n", ch.round(2).to_string(index=False))


# ---------------------------------------------------------------- 边际烟雾模型 (同 06 fe_ols: stratum 去均值 OLS, person 聚类)
def fe_ols(dd, y, xs, label):
    dd = dd.dropna(subset=[y] + xs)
    dd = dd[dd.groupby("stratum")[y].transform("size") >= 2]
    g = dd.groupby("stratum")
    Y = (dd[y] - g[y].transform("mean")).values
    X = np.column_stack([(dd[c] - g[c].transform("mean")).values for c in xs])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ Y)
    u = Y - X @ beta
    cp = pd.factorize(dd["person_id"])[0]
    G = cp.max() + 1; n, k = X.shape
    S = np.column_stack([np.bincount(cp, weights=X[:, j] * u, minlength=G) for j in range(k)])
    V = (G / (G - 1)) * ((n - 1) / (n - k)) * (XtX_inv @ (S.T @ S) @ XtX_inv)
    out = []
    for j, nm in enumerate(xs):
        se = np.sqrt(V[j, j])
        out.append({"model": label, "term": nm, "beta": beta[j], "se": se,
                    "lo": beta[j] - 1.96 * se, "hi": beta[j] + 1.96 * se,
                    "p": 2 * norm.sf(abs(beta[j] / se)), "n_nights": len(dd),
                    "n_persons": dd["person_id"].nunique()})
    return pd.DataFrame(out)


res = []
for dim, lab, m in LEVELS:
    o = fe_ols(d[m], "minute_asleep", ["sd", "tmean"], f"{dim} | {lab}")
    res.append(o)
    k = o[o.term == "sd"].iloc[0]
    log(f"{dim} | {lab}: beta={k.beta:+.3f} ({k.lo:+.2f}, {k.hi:+.2f}) p={k.p:.4f} "
        f"persons={k.n_persons:,} nights={k.n_nights:,}")
pd.concat(res).to_csv("rev5_subgroup_marginal.csv", index=False)

pd.DataFrame([{"key": "rows", "value": len(d)}, {"key": "persons", "value": d["person_id"].nunique()},
              {"key": "adi_p75_threshold_nights", "value": thr75},
              {"key": "runtime_s", "value": int(time.time() - T0)}]).to_csv("rev5_meta.csv", index=False)
log("ALL DONE — 导出前核对: (1) 6 个已有水平的 beta 应与 eTable 6 一致 (cardioresp -1.48, OSA -1.13, no OSA -0.33, "
    "<65 -0.65, >=65 -0.11, deprived -0.32); (2) rev5_*.csv 无 1-20 计数。")
