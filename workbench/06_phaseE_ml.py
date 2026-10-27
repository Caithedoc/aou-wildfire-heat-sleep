"""Phase E: ML 三件套 (2026-08-28)。在 AoU Jupyter 里跑, 读 VM 缓存, 不重下数据。

三件:
  E1 反事实预期睡眠基线 (梯度提升; 锚 Heffernan 2025 AJE ML-CITS)
      用非烟夜训练 "这个人这一晚该睡多久", 在烟夜看实测减预期的残差。
  E2 causal forest 脆弱性画像 (R-learner on 层内残差 + GATES 校准检验)
      X = 人级修饰因子 (年龄/慢病/剥夺指数/基线睡眠表型), 输出 tau(X) 分位组效应。
  E3 佩戴逆概率加权 (回应 Hurwitz 2026 佩戴偏倚框架)
      建模 "今晚会不会有有效睡眠记录", 用 stabilized IPW 重跑主效应。

纪律: 只导出聚合系数/分组效应表; 任何分组 n<20 不呈报。
"""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import statsmodels.api as sm
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.model_selection import KFold

DS = "wb-clever-rhubarb-4688.wildfire_sleep"
THRESH = 20
RNG = 20260828

# ---------------------------------------------------------------- 读缓存
# BigQuery 的 dbdate 扩展类型 pandas 读不回来, 逐列 to_numpy 绕开 (Phase F 实测)。
COLS = ["person_id", "sleep_date", "ym", "dow", "src_id", "zip3",
        "minute_asleep", "minute_awake", "sleep_efficiency",
        "smoke", "smoke_l1", "tmean"]
t = pq.read_table("analysis_nights.parquet")
have = [c for c in COLS if c in t.column_names]
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in have})
for c in ["minute_asleep", "minute_awake", "sleep_efficiency", "smoke", "smoke_l1", "tmean"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["date"] = pd.to_datetime(df["sleep_date"].astype(str), errors="coerce")
df["month"] = df["date"].dt.month
df["doy"] = df["date"].dt.dayofyear
df["year"] = df["date"].dt.year
df["smoke_day"] = (df["smoke"] > THRESH).astype(float)
df["stratum"] = (df["person_id"].astype(str) + "_" + df["ym"].astype(str)
                 + "_" + df["dow"].astype(str))
df = df.dropna(subset=["minute_asleep", "smoke", "tmean"])
print("rows:", len(df), "| persons:", df["person_id"].nunique())


def fe_ols(d, y, xs, label, weights=None):
    """层内去均值 OLS (可加权), 按 person 聚类 SE。"""
    d = d.dropna(subset=[y] + xs).copy()
    g = d.groupby("stratum")
    d = d[g[y].transform("size") >= 2]
    g = d.groupby("stratum")
    Y = (d[y] - g[y].transform("mean")).values
    X = np.column_stack([(d[c] - g[c].transform("mean")).values for c in xs])
    if weights is None:
        m = sm.OLS(Y, X).fit(cov_type="cluster",
                             cov_kwds={"groups": d["person_id"].values})
    else:
        w = d[weights].values
        m = sm.WLS(Y, X, weights=w).fit(cov_type="cluster",
                                        cov_kwds={"groups": d["person_id"].values})
    return pd.DataFrame({"model": label, "term": xs, "beta": m.params, "se": m.bse,
                         "lo": m.params - 1.96 * m.bse, "hi": m.params + 1.96 * m.bse,
                         "p": m.pvalues, "n_nights": int(m.nobs),
                         "n_persons": d["person_id"].nunique()})


results = []

# ================================================================ E1 反事实基线
# 训练集 = 非烟夜 (smoke <= 5, 干净对照), 特征 = 个体历史均值 + 季节 + 星期 + 温度。
print("\n=== E1 counterfactual baseline ===")
hist = (df[df["smoke"] <= 5].groupby("person_id")["minute_asleep"]
        .agg(["mean", "std", "count"]).rename(columns={
            "mean": "hist_mean", "std": "hist_sd", "count": "hist_n"}))
d1 = df.join(hist, on="person_id").dropna(subset=["hist_mean"])
d1 = d1[d1["hist_n"] >= 20]                      # 个体历史至少 20 夜才可外推
FEATS = ["hist_mean", "hist_sd", "tmean", "doy", "dow", "year"]
train = d1[d1["smoke"] <= 5]
mdl = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08,
                                    random_state=RNG)
mdl.fit(train[FEATS].values, train["minute_asleep"].values)
d1["pred"] = mdl.predict(d1[FEATS].values)
d1["resid"] = d1["minute_asleep"] - d1["pred"]
r2_in = 1 - ((train["minute_asleep"] - mdl.predict(train[FEATS].values)) ** 2).sum() / \
        ((train["minute_asleep"] - train["minute_asleep"].mean()) ** 2).sum()
print(f"baseline model in-sample R2 = {r2_in:.3f}")
# 烟夜残差 vs 非烟夜残差 (层内去均值, 与主分析同规格)
res_e1 = fe_ols(d1, "resid", ["smoke_day", "tmean"], "E1_residual_vs_expected")
print(res_e1.round(3).to_string(index=False))
results.append(res_e1)

# ================================================================ E2 causal forest
print("\n=== E2 heterogeneity (R-learner + GATES) ===")
cov = None
try:
    import subprocess
    from io import StringIO
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        "--max_rows=200000",
                        f"SELECT person_id, year_of_birth, deprivation_index, "
                        f"median_income, has_osa, has_insomnia, has_asthma, has_copd, "
                        f"has_t2dm, has_htn, has_ihd, has_hf FROM `{DS}.person_covars`"],
                       capture_output=True, text=True)
    cov = pd.read_csv(StringIO(r.stdout)) if r.returncode == 0 else None
except Exception as e:
    print("covars fetch failed:", e)

if cov is not None:
    d2 = df.merge(cov, on="person_id", how="inner")
    d2 = d2.join(hist, on="person_id")
    d2["age"] = d2["year"] - d2["year_of_birth"]
    MODS = ["age", "deprivation_index", "median_income", "hist_mean", "hist_sd",
            "has_osa", "has_insomnia", "has_asthma", "has_copd",
            "has_t2dm", "has_htn", "has_ihd", "has_hf"]
    d2 = d2.dropna(subset=["age", "hist_mean"]).copy()
    for m in MODS:
        d2[m] = pd.to_numeric(d2[m], errors="coerce")
    d2[MODS] = d2[MODS].fillna(d2[MODS].median())

    # 层内残差化 (Robinson): 去掉 person x 月 x 星期几 的一切固定效应
    g = d2.groupby("stratum")
    d2 = d2[g["minute_asleep"].transform("size") >= 2]
    g = d2.groupby("stratum")
    d2["y_t"] = d2["minute_asleep"] - g["minute_asleep"].transform("mean")
    d2["w_t"] = d2["smoke_day"] - g["smoke_day"].transform("mean")

    # 交叉拟合 tau(X): 用 (y_t * w_t) / var(w_t) 的伪结局做回归 (R-learner 权重形式)
    keep = d2["w_t"].abs() > 1e-8
    d2f = d2[keep].copy()
    d2f["pseudo"] = d2f["y_t"] / d2f["w_t"]
    d2f["wt"] = d2f["w_t"] ** 2
    tau_hat = np.zeros(len(d2f))
    persons = d2f["person_id"].values
    uniq = pd.unique(persons)
    kf = KFold(n_splits=3, shuffle=True, random_state=RNG)
    for tr_idx, te_idx in kf.split(uniq):                 # 按人切, 防同人跨折泄漏
        tr_p, te_p = set(uniq[tr_idx]), set(uniq[te_idx])
        m_tr = d2f["person_id"].isin(tr_p).values
        m_te = d2f["person_id"].isin(te_p).values
        f = HistGradientBoostingRegressor(max_iter=150, learning_rate=0.08,
                                          random_state=RNG)
        f.fit(d2f.loc[m_tr, MODS].values, d2f.loc[m_tr, "pseudo"].values,
              sample_weight=d2f.loc[m_tr, "wt"].values)
        tau_hat[m_te] = f.predict(d2f.loc[m_te, MODS].values)
    d2f["tau_hat"] = tau_hat

    # GATES: 按预测 tau 分五组, 每组用主分析同规格重估真实效应 (校准检验)
    d2f["tau_grp"] = pd.qcut(d2f["tau_hat"], 5, labels=False, duplicates="drop")
    gates = []
    for q, gq in d2f.groupby("tau_grp"):
        ids = set(gq["person_id"])
        sub = d2[d2["person_id"].isin(ids)]
        if sub["person_id"].nunique() < 20:
            continue
        rr = fe_ols(sub, "minute_asleep", ["smoke_day", "tmean"], f"E2_GATES_q{int(q)+1}")
        rr["pred_tau_mean"] = gq["tau_hat"].mean()
        gates.append(rr)
    if gates:
        gt = pd.concat(gates, ignore_index=True)
        print(gt[gt["term"] == "smoke_day"].round(3).to_string(index=False))
        results.append(gt)

    # 预设亚组 (不靠 ML, 直接报, 审稿人更认)
    SUBS = {"age_ge65": d2["age"] >= 65, "age_lt65": d2["age"] < 65,
            "osa": d2["has_osa"] == 1, "no_osa": d2["has_osa"] == 0,
            "cardioresp": (d2[["has_asthma", "has_copd", "has_ihd", "has_hf"]].sum(axis=1) > 0),
            "deprived_top25": d2["deprivation_index"] >= d2["deprivation_index"].quantile(0.75)}
    for name, mask in SUBS.items():
        sub = d2[mask]
        if sub["person_id"].nunique() < 20:
            continue
        results.append(fe_ols(sub, "minute_asleep", ["smoke_day", "tmean"],
                              "E2_subgroup_" + name))
    print(pd.concat([r for r in results if str(r["model"].iloc[0]).startswith("E2_subgroup")])
          .query("term == 'smoke_day'").round(3).to_string(index=False))

# ================================================================ E3 佩戴 IPW
print("\n=== E3 wear-day IPW ===")
try:
    import subprocess
    from io import StringIO
    # 每人 x 日: 是否有有效睡眠夜 (分析集内) vs 有佩戴但无睡眠 —— 用 wear_days 做全集
    q = f"""
    WITH span AS (
      SELECT person_id, MIN(sleep_date) d0, MAX(sleep_date) d1
      FROM `{DS}.cohort_nights` GROUP BY person_id
    ),
    w AS (
      SELECT w.person_id, w.wear_date AS date
      FROM `{DS}.wear_days` w JOIN span s USING (person_id)
      WHERE w.wear_date BETWEEN s.d0 AND s.d1
    )
    SELECT w.person_id, w.date,
           IF(c.person_id IS NULL, 0, 1) AS has_sleep,
           e.smokePM_popwt AS smoke, e.tmean_popwt AS tmean
    FROM w
    LEFT JOIN `{DS}.cohort_nights` c
      ON c.person_id = w.person_id AND c.sleep_date = w.date
    JOIN `{DS}.person_zip3` z ON z.person_id = w.person_id
    JOIN `{DS}.exposure_zip3_daily` e ON e.zip3 = z.zip3 AND e.date = w.date
    """
    r = subprocess.run(["bq", "query", "--nouse_legacy_sql", "--format=csv",
                        "--max_rows=100000000", q], capture_output=True, text=True)
    wd = pd.read_csv(StringIO(r.stdout)) if r.returncode == 0 else None
except Exception as e:
    wd = None
    print("wear panel fetch failed:", e)

if wd is not None and len(wd):
    wd["date"] = pd.to_datetime(wd["date"])
    wd["smoke_day"] = (pd.to_numeric(wd["smoke"], errors="coerce") > THRESH).astype(float)
    wd["doy"] = wd["date"].dt.dayofyear
    wd["dow"] = wd["date"].dt.dayofweek
    wd["year"] = wd["date"].dt.year
    wd["tmean"] = pd.to_numeric(wd["tmean"], errors="coerce")
    wd = wd.dropna(subset=["tmean"])
    print("wear-day panel rows:", len(wd),
          "| has_sleep rate:", round(wd["has_sleep"].mean(), 3))
    # 第二结局: 烟日是否降低"有睡眠记录"的概率 (缺失机制本身)
    ph = wd.groupby("person_id")["has_sleep"].mean().rename("p_hist")
    wd = wd.join(ph, on="person_id")
    F2 = ["p_hist", "tmean", "doy", "dow", "year", "smoke_day"]
    clf = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.08,
                                         random_state=RNG)
    clf.fit(wd[F2].values, wd["has_sleep"].values)
    wd["p_obs"] = clf.predict_proba(wd[F2].values)[:, 1].clip(0.05, 0.999)
    wd["sw"] = wd["has_sleep"].mean() / wd["p_obs"]          # stabilized IPW
    # 烟日 vs 非烟日的佩戴/记录差 (描述性, 直接可报)
    obs_rate = wd.groupby("smoke_day")["has_sleep"].agg(["mean", "size"])
    print("has_sleep rate by smoke day:\n", obs_rate.round(4).to_string())
    # 把权重挂回分析集重跑主效应
    wsub = wd.loc[wd["has_sleep"] == 1, ["person_id", "date", "sw"]]
    d3 = df.merge(wsub, left_on=["person_id", "date"], right_on=["person_id", "date"],
                  how="inner")
    print("weighted analysis rows:", len(d3))
    results.append(fe_ols(d3, "minute_asleep", ["smoke_day", "tmean"],
                          "E3_IPW_weighted", weights="sw"))
    results.append(fe_ols(d3, "minute_asleep", ["smoke_day", "tmean"],
                          "E3_unweighted_same_sample"))

# ================================================================ E4 烟 x 热 联合效应
# 不用 WQS/BKMR (为人际间化学混合物设计, 与个体内急性设计不兼容, 且 index p 值种子不稳)。
# 真问题 = 烟与热同季共现时的联合/交互效应, 在 case-crossover 框架内直接可答。
print("\n=== E4 smoke x heat joint effect ===")
d4 = df.copy()
# 热日 = 该 zip3 该日历月温度分布的 90 分位以上 (个体内比较, 相对定义)
d4["zm"] = d4["zip3"].astype(str) + "_" + d4["month"].astype(str)
p90 = d4.groupby("zm")["tmean"].transform(lambda s: s.quantile(0.90))
d4["heat_day"] = (d4["tmean"] >= p90).astype(float)
d4["smoke_x_heat"] = d4["smoke_day"] * d4["heat_day"]
print("heat day share:", round(d4["heat_day"].mean(), 3),
      "| both share:", round(d4["smoke_x_heat"].mean(), 4))

# (a) 乘法/加法尺度交互 (连续结局下即加法尺度)
results.append(fe_ols(d4, "minute_asleep",
                      ["smoke_day", "heat_day", "smoke_x_heat", "tmean"],
                      "E4_interaction"))
# (b) 四格联合暴露 (参照 = 都无), 最直观, 可直接做表
d4["g_heat_only"] = ((d4["heat_day"] == 1) & (d4["smoke_day"] == 0)).astype(float)
d4["g_smoke_only"] = ((d4["heat_day"] == 0) & (d4["smoke_day"] == 1)).astype(float)
d4["g_both"] = ((d4["heat_day"] == 1) & (d4["smoke_day"] == 1)).astype(float)
n_both = int(d4["g_both"].sum())
print("nights in both-exposure cell:", n_both)
if n_both >= 20:
    results.append(fe_ols(d4, "minute_asleep",
                          ["g_heat_only", "g_smoke_only", "g_both"],
                          "E4_joint_4cell"))

out = pd.concat(results, ignore_index=True)
out.to_csv("phaseE_results.csv", index=False)
print("\nsaved phaseE_results.csv")
print(out[out["term"] == "smoke_day"].round(3).to_string(index=False))
