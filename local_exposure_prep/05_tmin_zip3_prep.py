"""Phase C 补充: PRISM tmin -> zip3-日 面板 (2026-08-30, 本地跑)。

与 03_smokepm_zip3_prep.py stage_temp 完全同构 (同一批 ZCTA 质心、同一人口权重、
同一 CONUS 掩膜), 只把变量换成日最低温。产出:
    out/tmin_zip3_daily.parquet
    out/tmin_zip3_daily.csv.gz   <- 上传到 AoU workspace bucket, 供 12_recluster_tempflex.py 合并

跑法:
    .venv/bin/python 05_tmin_zip3_prep.py
前置: download_prism_tmin.sh 已跑完 (raw/prism_tmin/ 应有 2,925 个 zip)。
"""

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
RAW = HERE / "raw"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)


def main():
    import rasterio

    pop = pd.read_csv(RAW / "zcta_population.csv", dtype={"zcta": str})
    pop["zip3"] = pop["zcta"].str[:3]
    gaz = pd.read_csv(RAW / "2019_Gaz_zcta_national.txt", sep="\t", dtype={"GEOID": str})
    gaz.columns = [c.strip() for c in gaz.columns]
    gaz = gaz.rename(columns={"GEOID": "zcta", "INTPTLAT": "lat", "INTPTLONG": "lon"})
    pts = gaz[["zcta", "lat", "lon"]].merge(pop[["zcta", "pop", "zip3"]], on="zcta", how="inner")
    pts = pts[pts["pop"] > 0].reset_index(drop=True)
    print(f"采样点 {len(pts):,} 个 ZCTA 质心")

    pdir = RAW / "prism_tmin"
    zips = sorted(pdir.glob("prism_tmin_us_25m_*.zip"))
    print(f"PRISM tmin 文件 {len(zips):,} 个")
    assert len(zips) >= 2900, "文件数不足, 先确认 download_prism_tmin.sh 跑完且 fail.log 为空"

    rows_idx = cols_idx = None
    recs = []
    for i, zp in enumerate(zips):
        ymd = zp.stem.split("_")[-1]
        tif = f"zip://{zp}!{zp.stem}.tif"
        try:
            with rasterio.open(tif) as src:
                if rows_idx is None:
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
        w = pts["pop"].values * ok
        num = pd.Series(t * w).groupby(pts["zip3"]).sum()
        den = pd.Series(w).groupby(pts["zip3"]).sum()
        tm = (num / den.replace(0, np.nan)).dropna()
        recs.append(pd.DataFrame({"zip3": tm.index, "date": pd.Timestamp(ymd),
                                  "tmin_popwt": tm.values}))
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(zips)}")
    temp = pd.concat(recs, ignore_index=True)
    temp.to_parquet(OUT / "tmin_zip3_daily.parquet", index=False)
    csv = temp.copy()
    csv["date"] = csv["date"].dt.strftime("%Y-%m-%d")
    csv.to_csv(OUT / "tmin_zip3_daily.csv.gz", index=False, compression="gzip")
    print(f"tmin 面板 {len(temp):,} 行, zip3 数 {temp['zip3'].nunique()}")

    # 一致性抽检: 与 tmean 面板同 zip3-日, tmin 应恒 <= tmean (数值容差 0.5C)
    try:
        tmean = pd.read_parquet(OUT / "tmean_zip3_daily.parquet")
        m = temp.merge(tmean, on=["zip3", "date"], how="inner")
        bad = (m["tmin_popwt"] > m["tmean_popwt"] + 0.5).mean()
        print(f"抽检: 合并 {len(m):,} 行; tmin>tmean+0.5C 比例 {bad:.6f} (应≈0)")
        print("tmin 均值 {:.2f}C, tmean 均值 {:.2f}C, 平均日较差/2 ≈ {:.2f}C".format(
            m["tmin_popwt"].mean(), m["tmean_popwt"].mean(),
            (m["tmean_popwt"] - m["tmin_popwt"]).mean()))
    except Exception as e:
        print("tmean 对照抽检失败:", e)


if __name__ == "__main__":
    main()
