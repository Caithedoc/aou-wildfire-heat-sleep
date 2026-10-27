"""07: ZCTA -> ZIP3 聚合误差量化 (2026-09-22, 本地跑, 全部公开数据, 不需要 AoU)。

背景: AoU 只给 3 位 ZIP, 暴露按 2020 人口加权从 ZCTA 聚到 ZIP3-日。审稿人最可能追问
ZIP3 内空间异质性造成的暴露误分类。本脚本回答:
  (1) 以 ZCTA 为"真值"、ZIP3 为"赋值", 人口加权的人-日层面一致性
      (灵敏度 / PPV / 特异度 / kappa), 分烟雾日 / 高温日 / 复合日, 分地区与年份;
  (2) 每个 ZIP3 的烟雾分类异质性指数 H, 供 VM 上做"限低异质性 ZIP3"敏感性分析。

事前规定 (写于看到任何 AoU 结果之前, 2026-09-22):
  H_smoke(z) = 误分类人-日 / (ZCTA 或 ZIP3 任一判为烟雾日的人-日)   [= 1 - 人口加权 Jaccard]
  限定分析 R1: H_smoke <= 分析 ZIP3 的人口加权中位数 (在 VM 上改用分析夜加权中位数, 同时报两者)
  限定分析 R2: H_smoke 最低三分位
  限定分析 R3: 剔除 H_smoke 最高三分位
  从未出现烟雾日的 ZIP3 (分母为 0) H 记 0, 但不影响交互项 (无烟雾夜)。
  估计量与主分析一致: 交互参数化 sd + hd + sxh + tmean 样条 (M3), 四种聚类 + ZIP3 WCB。

分段:
    .venv/bin/python 07_zip3_aggregation_error.py ztemp    # PRISM tmean 在每个 ZCTA 质心取值 -> 矩阵
    .venv/bin/python 07_zip3_aggregation_error.py zsmoke   # rds -> ZCTA x 日 烟雾矩阵
    .venv/bin/python 07_zip3_aggregation_error.py analyze  # 一致性指标 + ZIP3 异质性表
输出 out/agg_error/。

口径说明: 分析里 P90 是 zip3 x 月 在"分析夜"上取分位 (夜加权), 这里用 2016-2023 全部日取分位
(日加权), 两者差别仅在权重, 已在输出 meta 注明。
"""

import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
RAW = HERE / "raw"
OUT = HERE / "out"
AE = OUT / "agg_error"
AE.mkdir(parents=True, exist_ok=True)

WIN_START, WIN_END = "2016-01-01", "2023-12-31"
DAYS = pd.date_range(WIN_START, WIN_END, freq="D")
SMOKE_THRS = (10, 20, 35)
MAIN_THR = 20
HEAT_Q = 0.90


def region_of(zip3: pd.Series) -> pd.Series:
    # 同 workbench/11_revision.py:308 — ZIP 首位 0-2 East, 3-7 Central, 8-9 West
    f = zip3.str[0].astype(int)
    return pd.Series(np.select([f <= 2, f <= 7], ["East", "Central"], "West"), index=zip3.index)


def universe() -> pd.DataFrame:
    """ZCTA 宇宙: 2020 人口 > 0 且 ZIP3 在 PRISM 覆盖的 CONUS 分析面板里 (880 个)。"""
    pop = pd.read_csv(RAW / "zcta_population.csv", dtype={"zcta": str})
    pop["zip3"] = pop["zcta"].str[:3]
    conus = set(pd.read_parquet(OUT / "tmean_zip3_daily.parquet", columns=["zip3"])["zip3"].unique())
    u = pop[(pop["pop"] > 0) & pop["zip3"].isin(conus)].sort_values("zcta").reset_index(drop=True)
    return u


# ---------------------------------------------------------------- ZCTA 温度矩阵
def stage_ztemp():
    import rasterio
    u = universe()
    gaz = pd.read_csv(RAW / "2019_Gaz_zcta_national.txt", sep="\t", dtype={"GEOID": str})
    gaz.columns = [c.strip() for c in gaz.columns]
    gaz = gaz.rename(columns={"GEOID": "zcta", "INTPTLAT": "lat", "INTPTLONG": "lon"})
    pts = u.merge(gaz[["zcta", "lat", "lon"]], on="zcta", how="left")
    has = pts["lat"].notna().values
    print(f"ZCTA {len(pts):,}; 有质心 {has.sum():,} (人口占比 "
          f"{pts.loc[has, 'pop'].sum() / pts['pop'].sum():.4f})")
    M = np.full((len(pts), len(DAYS)), np.nan, dtype=np.float32)
    ri = ci = None
    for j, day in enumerate(DAYS):
        ymd = day.strftime("%Y%m%d")
        zp = RAW / "prism_tmean" / f"prism_tmean_us_25m_{ymd}.zip"
        if not zp.exists():
            print("MISSING", ymd); continue
        with rasterio.open(f"zip://{zp}!{zp.stem}.tif") as src:
            if ri is None:
                tr = src.transform
                ci = np.clip(np.floor((pts["lon"].fillna(0).values - tr.c) / tr.a).astype(int), 0, src.width - 1)
                ri = np.clip(np.floor((pts["lat"].fillna(0).values - tr.f) / tr.e).astype(int), 0, src.height - 1)
            arr = src.read(1)
            nod = src.nodata if src.nodata is not None else -9999.0
        t = arr[ri, ci].astype(np.float32)
        t[(t == nod) | ~has] = np.nan
        M[:, j] = t
        if (j + 1) % 365 == 0:
            print(f"  {j+1}/{len(DAYS)}", flush=True)
    np.save(OUT / "zcta_tmean_matrix.npy", M)
    pts[["zcta", "zip3", "pop"]].to_csv(OUT / "zcta_matrix_index.csv", index=False)
    print("saved zcta_tmean_matrix.npy", M.shape, "NaN 比例", float(np.isnan(M).mean()))


# ---------------------------------------------------------------- ZCTA 烟雾矩阵
def stage_zsmoke():
    import pyreadr
    u = universe()
    df = next(iter(pyreadr.read_r(str(RAW / "smokePM2pt5_predictions_daily_zcta_20060101-20231231.rds")).values()))
    cols = {c.lower(): c for c in df.columns}
    zcol = cols.get("zcta") or cols.get("zip") or cols.get("geoid")
    df = df[[zcol, cols["date"], cols.get("smokepm_pred") or cols.get("smokepm")]]
    df.columns = ["zcta", "date", "smoke"]
    df["zcta"] = df["zcta"].astype(str).str.zfill(5)
    df["date"] = pd.to_datetime(df["date"])
    df = df[(df["date"] >= WIN_START) & (df["date"] <= WIN_END)]
    zi = pd.Series(np.arange(len(u)), index=u["zcta"])
    df = df[df["zcta"].isin(zi.index)]
    S = np.zeros((len(u), len(DAYS)), dtype=np.float32)   # 产品语义: 无烟羽 = 0
    S[zi[df["zcta"]].values, (df["date"] - DAYS[0]).dt.days.values] = df["smoke"].values
    np.save(OUT / "zcta_smoke_matrix.npy", S)
    u[["zcta", "zip3", "pop"]].to_csv(OUT / "zcta_smoke_index.csv", index=False)
    print("saved zcta_smoke_matrix.npy", S.shape, "非零比例", float((S > 0).mean()))


# ---------------------------------------------------------------- 分析
def wsum_by(codes, x, G):
    return np.bincount(codes, weights=x, minlength=G)


def agreement(tru, asg, w):
    """人口加权人-日 2x2: tru = ZCTA 判定, asg = ZIP3 赋值 (bool 矩阵), w 与之同形或可广播。"""
    tp = float((w * (tru & asg)).sum()); fp = float((w * (~tru & asg)).sum())
    fn = float((w * (tru & ~asg)).sum()); tn = float((w * (~tru & ~asg)).sum())
    n = tp + fp + fn + tn
    po = (tp + tn) / n
    pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n ** 2
    return {"prev_zcta": (tp + fn) / n, "prev_zip3": (tp + fp) / n,
            "sensitivity": tp / (tp + fn) if tp + fn else np.nan,
            "ppv": tp / (tp + fp) if tp + fp else np.nan,
            "specificity": tn / (tn + fp) if tn + fp else np.nan,
            "kappa": (po - pe) / (1 - pe) if pe < 1 else np.nan,
            "misclassified_share_of_either": (fp + fn) / (tp + fp + fn) if tp + fp + fn else np.nan,
            "person_days_either_M": (tp + fp + fn) / 1e6}


def monthly_p90(X, months):
    """X: (units, days); 每行每个日历月的 P90 (nan 忽略) -> 同形阈值矩阵。"""
    thr = np.full(X.shape, np.nan, dtype=np.float32)
    for m in range(1, 13):
        cols = np.where(months == m)[0]
        q = np.nanquantile(X[:, cols], HEAT_Q, axis=1)
        thr[:, cols] = q[:, None]
    return thr


def stage_analyze():
    S = np.load(OUT / "zcta_smoke_matrix.npy")
    idx = pd.read_csv(OUT / "zcta_smoke_index.csv", dtype={"zcta": str, "zip3": str})
    T = np.load(OUT / "zcta_tmean_matrix.npy")
    tidx = pd.read_csv(OUT / "zcta_matrix_index.csv", dtype={"zcta": str, "zip3": str})
    assert (idx["zcta"].values == tidx["zcta"].values).all(), "两矩阵 ZCTA 顺序不一致"
    pop = idx["pop"].values.astype(np.float64)
    z3, z3_codes = np.unique(idx["zip3"].values, return_inverse=True)
    G = len(z3)
    months = DAYS.month.values
    years = DAYS.year.values

    # ---- ZIP3 烟雾: 按 03 同口径重算 (分母 = ZIP3 全部有人口 ZCTA) 并与已用面板核对
    zpop = wsum_by(z3_codes, pop, G)
    Z = np.vstack([wsum_by(z3_codes, S[:, j] * pop, G) for j in range(S.shape[1])]).T / zpop[:, None]
    panel = pd.read_parquet(OUT / "smoke_zip3_daily.parquet")
    panel = panel[panel["zip3"].isin(z3)]
    pv = panel.pivot(index="zip3", columns="date", values="smokePM_popwt").reindex(index=z3, columns=DAYS)
    diff = np.nanmax(np.abs(pv.values - Z))
    print(f"ZIP3 烟雾重算 vs 分析面板 最大绝对差 = {diff:.2e} µg/m³")
    assert diff < 0.05, "重算与分析面板不一致, 停止"

    # ---- ZIP3 温度 = 分析面板值; ZCTA 温度 = 质心值
    tp = pd.read_parquet(OUT / "tmean_zip3_daily.parquet")
    TZ = tp.pivot(index="zip3", columns="date", values="tmean_popwt").reindex(index=z3, columns=DAYS).values
    t_ok = ~np.isnan(T).all(axis=1)
    print(f"有温度的 ZCTA 人口占比 {pop[t_ok].sum() / pop.sum():.4f}")
    heat_z3 = TZ >= monthly_p90(TZ, months)
    heat_zc = T >= monthly_p90(T, months)

    W = pop[:, None]
    rows, zrows = [], []
    reg = region_of(pd.Series(z3)).values
    reg_zc = reg[z3_codes]

    def emit(label, tru, asg, mask_rows=None, mask_cols=None, extra=None):
        tr_, as_, w_ = tru, asg, np.broadcast_to(W, tru.shape)
        if mask_rows is not None:
            tr_, as_, w_ = tr_[mask_rows], as_[mask_rows], w_[mask_rows]
        if mask_cols is not None:
            tr_, as_, w_ = tr_[:, mask_cols], as_[:, mask_cols], w_[:, mask_cols]
        r = {"exposure": label, **(extra or {}), **agreement(tr_, as_, w_)}
        rows.append(r)

    for thr in SMOKE_THRS:
        tru = S > thr
        asg = (Z > thr)[z3_codes]
        emit(f"smoke>{thr}", tru, asg, extra={"stratum": "all"})
        if thr == MAIN_THR:
            for rg in ("East", "Central", "West"):
                emit(f"smoke>{thr}", tru, asg, mask_rows=reg_zc == rg, extra={"stratum": rg})
            for y in range(2016, 2024):
                emit(f"smoke>{thr}", tru, asg, mask_cols=years == y, extra={"stratum": str(y)})

    ok = t_ok
    heat_asg = heat_z3[z3_codes]
    emit("heat_P90", heat_zc, heat_asg, mask_rows=ok, extra={"stratum": "all (ZCTA with centroid)"})
    comp_tru = (S > MAIN_THR) & heat_zc
    comp_asg = ((Z > MAIN_THR) & heat_z3)[z3_codes]
    emit("compound", comp_tru, comp_asg, mask_rows=ok, extra={"stratum": "all (ZCTA with centroid)"})
    for rg in ("East", "Central", "West"):
        m = ok & (reg_zc == rg)
        emit("compound", comp_tru, comp_asg, mask_rows=m, extra={"stratum": rg})
        emit("heat_P90", heat_zc, heat_asg, mask_rows=m, extra={"stratum": rg})

    ag = pd.DataFrame(rows)
    ag.to_csv(AE / "agreement_summary.csv", index=False)
    print(ag.round(4).to_string(index=False))

    # ---- 连续尺度: ZIP3 烟雾日上 ZCTA 值的人口加权组内 SD 与 ZCTA-ZIP3 相关
    smoke_days = (Z > MAIN_THR)[z3_codes]
    Zc = Z[z3_codes]
    dev = (S - Zc) ** 2
    wsd = np.sqrt(float((W * dev * smoke_days).sum()) / float((W * smoke_days).sum()))
    mask_any = (S > 0) | (Zc > 0)
    x, y, w = S[mask_any], Zc[mask_any], np.broadcast_to(W, S.shape)[mask_any]
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    r = np.sum(w * (x - mx) * (y - my)) / np.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2))
    cont = pd.DataFrame([{"within_zip3_sd_on_zip3_smoke_days": wsd,
                          "pop_weighted_r_zcta_vs_zip3_any_smoke": r,
                          "zcta_n": len(idx), "zip3_n": G,
                          "zip3_recompute_max_abs_diff": diff,
                          "p90_weighting": "days 2016-2023 (analysis uses nights)"}])
    cont.to_csv(AE / "continuous_meta.csv", index=False)
    print(cont.round(4).to_string(index=False))

    # ---- 每个 ZIP3 的异质性指数
    tru, asg = S > MAIN_THR, (Z > MAIN_THR)[z3_codes]
    mis = wsum_by(z3_codes, (W * (tru ^ asg)).sum(axis=1), G)
    eit = wsum_by(z3_codes, (W * (tru | asg)).sum(axis=1), G)
    cmis = wsum_by(z3_codes[ok], (W[ok] * (comp_tru[ok] ^ comp_asg[ok])).sum(axis=1), G)
    ceit = wsum_by(z3_codes[ok], (W[ok] * (comp_tru[ok] | comp_asg[ok])).sum(axis=1), G)
    nz = np.bincount(z3_codes, minlength=G)
    top = pd.Series(pop).groupby(z3_codes).max().values / zpop
    zt = pd.DataFrame({"zip3": z3, "region": reg, "n_zcta": nz, "pop2020": zpop.astype(int),
                       "largest_zcta_pop_share": top,
                       "zip3_smoke_days": (Z > MAIN_THR).sum(axis=1),
                       "zip3_compound_days": ((Z > MAIN_THR) & heat_z3).sum(axis=1),
                       "H_smoke": np.where(eit > 0, mis / np.maximum(eit, 1), 0.0),
                       "H_compound": np.where(ceit > 0, cmis / np.maximum(ceit, 1), 0.0)})
    ever = zt["zip3_smoke_days"] > 0
    zs = zt[ever].sort_values("H_smoke")
    cw = zs["pop2020"].cumsum() / zs["pop2020"].sum()
    med = float(zs["H_smoke"].values[np.searchsorted(cw.values, 0.5)])
    t1, t2 = np.quantile(zs["H_smoke"], [1 / 3, 2 / 3])
    zt["R1_le_popwt_median"] = (zt["H_smoke"] <= med).astype(int)
    zt["R2_bottom_tercile"] = (zt["H_smoke"] <= t1).astype(int)
    zt["R3_drop_top_tercile"] = (zt["H_smoke"] <= t2).astype(int)
    zt.to_csv(AE / "zip3_heterogeneity.csv", index=False)
    print(f"\nH_smoke (有烟雾日的 {ever.sum()} 个 ZIP3): 人口加权中位数 {med:.3f}; "
          f"三分位 {t1:.3f} / {t2:.3f}")
    print(zt[ever].groupby("region")["H_smoke"].describe().round(3).to_string())
    print("\n单 ZCTA 的 ZIP3 数:", int((zt["n_zcta"] == 1).sum()))
    print("H_smoke vs n_zcta 相关:", round(zt.loc[ever, ["H_smoke", "n_zcta"]].corr().iloc[0, 1], 3))


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("ztemp", "all"):
        stage_ztemp()
    if stage in ("zsmoke", "all"):
        stage_zsmoke()
    if stage in ("analyze", "all"):
        stage_analyze()
