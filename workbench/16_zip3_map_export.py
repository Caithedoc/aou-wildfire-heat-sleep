"""16_zip3_map_export: 按 3-digit ZIP 汇总分析样本, 供 Fig (Lu 2025 Nat Health Fig 1 风格气泡地图) 使用 (2026-09-28).

在 AoU Jupyter VM 里跑: `python -u 16_zip3_map_export.py`  (读同目录缓存 analysis_nights.parquet, 秒级)
输出 rev6_zip3_map.csv (聚合; 1-20 人数/夜数置 NaN, 且该行其余统计量一并置 NaN):
  zip3, n_persons, n_nights, tst_mean, tst_sd, smoke_nights, smoke_pct, tmean_mean
"""
import numpy as np, pandas as pd, pyarrow.parquet as pq
MIN_CELL = 20
COLS = ["person_id", "zip3", "minute_asleep", "smoke", "tmean"]
t = pq.read_table("analysis_nights.parquet")
df = pd.DataFrame({c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS})
for c in ["minute_asleep", "smoke", "tmean"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
df["zip3"] = df["zip3"].astype(str).str.zfill(3)
g = df.groupby("zip3")
out = pd.DataFrame({
    "n_persons": g["person_id"].nunique(),
    "n_nights": g.size(),
    "tst_mean": g["minute_asleep"].mean(),
    "tst_sd": g["minute_asleep"].std(),
    "smoke_nights": g["smoke"].sum(),
    "tmean_mean": g["tmean"].mean(),
}).reset_index()
out["smoke_pct"] = 100 * out["smoke_nights"] / out["n_nights"]
small = (out["n_persons"] < MIN_CELL) | (out["n_nights"] < MIN_CELL) | ((out["smoke_nights"] > 0) & (out["smoke_nights"] < MIN_CELL))
out.loc[small, out.columns.difference(["zip3"])] = np.nan
out.to_csv("rev6_zip3_map.csv", index=False)
print(out.describe()); print("suppressed rows:", int(small.sum()), "of", len(out))
