"""Phase C: zip3-日 暴露面板制备 (2026-08-28 完整实现, 本地跑)。

三段式, 命令行参数选段:
    .venv/bin/python 03_smokepm_zip3_prep.py smoke   # 烟 PM2.5 -> zip3-日 (人口加权, 补零)
    .venv/bin/python 03_smokepm_zip3_prep.py temp    # PRISM tmean -> zip3-日 (质心采样, 人口加权)
    .venv/bin/python 03_smokepm_zip3_prep.py merge   # 合并 -> out/smoke_temp_zip3_daily.csv

数据源 (raw/):
- smokePM2pt5_predictions_daily_zcta_20060101-20231231.rds (Stanford ECHO Lab v2.0,
  Childs 2022 EST 10.1021/acs.est.2c02934; 只在有烟羽的 zcta-日 给非零值)
- zcta_population.csv (2020 Census DHC P1_001N, 33,774 ZCTA)
- 2019_Gaz_zcta_national.txt (Census gazetteer, ZCTA 质心经纬度)
- prism_tmean/prism_tmean_us_25m_YYYYMMDD.zip (PRISM an 4km GeoTIFF, download_prism.sh)

面板范围: 2016-01-01 至 2023-12-31; zip3 集合 = PRISM 有值的 CONUS zip3
(AK/HI/PR 无 PRISM 覆盖, 合并时自动剔除, 分析随之限 CONUS)。
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
RAW = HERE / "raw"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

WIN_START, WIN_END = "2016-01-01", "2023-12-31"
DAYS = pd.date_range(WIN_START, WIN_END, freq="D")


def load_pop() -> pd.DataFrame:
    pop = pd.read_csv(RAW / "zcta_population.csv", dtype={"zcta": str})
    pop["zip3"] = pop["zcta"].str[:3]
    return pop


# ---------------------------------------------------------------- smoke
def stage_smoke():
    import pyreadr
    f = RAW / "smokePM2pt5_predictions_daily_zcta_20060101-20231231.rds"
    print("reading rds (266MB, 需要几分钟) ...")
    df = next(iter(pyreadr.read_r(str(f)).values()))
    # 列名以实际为准: 预期 zcta/zip + date + smokePM_pred
    cols = {c.lower(): c for c in df.columns}
    zcol = cols.get("zcta") or cols.get("zip") or cols.get("geoid")
    vcol = cols.get("smokepm_pred") or cols.get("smokepm")
    dcol = cols.get("date")
    assert zcol and vcol and dcol, f"未识别列名: {list(df.columns)}"
    df = df[[zcol, dcol, vcol]].rename(columns={zcol: "zcta", dcol: "date", vcol: "smoke"})
    df["zcta"] = df["zcta"].astype(str).str.zfill(5)
    df["date"] = pd.to_datetime(df["date"])
    print("rds 覆盖:", df["date"].min(), "->", df["date"].max(), f"; 行数 {len(df):,}")
    df = df[(df["date"] >= WIN_START) & (df["date"] <= WIN_END)]
    print("窗内烟行数:", f"{len(df):,}")

    pop = load_pop()
    zip3_pop = pop.groupby("zip3")["pop"].sum().rename("zip3_pop")
    m = df.merge(pop[["zcta", "pop", "zip3"]], on="zcta", how="inner")
    g = (m.assign(w=m["smoke"] * m["pop"])
           .groupby(["zip3", "date"], as_index=False)[["w"]].sum()
           .join(zip3_pop, on="zip3"))
    g["smokePM_popwt"] = g["w"] / g["zip3_pop"].clip(lower=1)
    g = g[["zip3", "date", "smokePM_popwt"]]
    # 补零: 全 zip3 x 全日期 (无烟羽 = 0, 这是该产品的语义)
    all_zip3 = sorted(zip3_pop.index)
    idx = pd.MultiIndex.from_product([all_zip3, DAYS], names=["zip3", "date"])
    panel = (g.set_index(["zip3", "date"]).reindex(idx, fill_value=0.0)
               .reset_index())
    panel.to_parquet(OUT / "smoke_zip3_daily.parquet", index=False)
    n_event = int((panel["smokePM_popwt"] > 20).sum())
    print(f"smoke 面板 {len(panel):,} 行; 烟日(>20) {n_event:,} zip3-日; "
          f"烟日(>10) {int((panel['smokePM_popwt'] > 10).sum()):,}")


# ---------------------------------------------------------------- temp
def stage_temp():
    import zipfile

    import rasterio
    pop = load_pop()
    gaz = pd.read_csv(RAW / "2019_Gaz_zcta_national.txt", sep="\t", dtype={"GEOID": str})
    gaz.columns = [c.strip() for c in gaz.columns]
    gaz = gaz.rename(columns={"GEOID": "zcta", "INTPTLAT": "lat", "INTPTLONG": "lon"})
    pts = gaz[["zcta", "lat", "lon"]].merge(pop[["zcta", "pop", "zip3"]], on="zcta", how="inner")
    pts = pts[pts["pop"] > 0].reset_index(drop=True)
    print(f"采样点 {len(pts):,} 个 ZCTA 质心")

    pdir = RAW / "prism_tmean"
    zips = sorted(pdir.glob("prism_tmean_us_25m_*.zip"))
    print(f"PRISM 文件 {len(zips):,} 个")

    rows_idx = cols_idx = None
    valid = None
    recs = []
    for i, zp in enumerate(zips):
        ymd = zp.stem.split("_")[-1]
        tif = f"zip://{zp}!{zp.stem}.tif"
        try:
            with rasterio.open(tif) as src:
                if rows_idx is None:
                    # rasterio 1.5 的 index() 不再接受数组, 用仿射变换手算
                    tr = src.transform
                    cols_idx = np.clip(np.floor((pts["lon"].values - tr.c) / tr.a).astype(int),
                                       0, src.width - 1)
                    rows_idx = np.clip(np.floor((pts["lat"].values - tr.f) / tr.e).astype(int),
                                       0, src.height - 1)
                arr = src.read(1)
                nod = src.nodata if src.nodata is not None else -9999.0
        except Exception as e:
            print("SKIP", ymd, e)
            continue
        t = arr[rows_idx, cols_idx].astype(float)
        ok = t != nod
        if valid is None:
            valid = ok  # 网格恒定, 有效掩膜首日锁定 (CONUS 内质心)
        w = pts["pop"].values * ok
        num = pd.Series(t * w).groupby(pts["zip3"]).sum()
        den = pd.Series(w).groupby(pts["zip3"]).sum()
        tm = (num / den.replace(0, np.nan)).dropna()
        recs.append(pd.DataFrame({"zip3": tm.index, "date": pd.Timestamp(ymd),
                                  "tmean_popwt": tm.values}))
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(zips)}")
    temp = pd.concat(recs, ignore_index=True)
    temp.to_parquet(OUT / "tmean_zip3_daily.parquet", index=False)
    print(f"temp 面板 {len(temp):,} 行, zip3 数 {temp['zip3'].nunique()}")


# ---------------------------------------------------------------- merge
def stage_merge():
    smoke = pd.read_parquet(OUT / "smoke_zip3_daily.parquet")
    temp = pd.read_parquet(OUT / "tmean_zip3_daily.parquet")
    panel = temp.merge(smoke, on=["zip3", "date"], how="inner")  # 限 PRISM 覆盖 zip3
    panel = panel.sort_values(["zip3", "date"])
    panel["date"] = panel["date"].dt.strftime("%Y-%m-%d")
    f = OUT / "smoke_temp_zip3_daily.csv"
    panel[["zip3", "date", "smokePM_popwt", "tmean_popwt"]].to_csv(f, index=False)
    print(f"最终面板 {len(panel):,} 行 -> {f}")
    print("烟日(>20):", int((panel["smokePM_popwt"] > 20).sum()),
          "; zip3:", panel["zip3"].nunique())


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("smoke", "all"):
        stage_smoke()
    if stage in ("temp", "all"):
        stage_temp()
    if stage in ("merge", "all"):
        stage_merge()
