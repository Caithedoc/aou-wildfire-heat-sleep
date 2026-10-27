"""Bubble map of participants by 3-digit ZIP, styled after Lu et al. 2025 Nat Health Fig 1
(state outlines, bubble size = participants, colour = mean total sleep time, viridis, size legend).
Run: /usr/bin/python3 figures/make_fig_map_zip3.py [rev6_zip3_map.csv]
Without the VM export it runs in PROXY mode on local exposure data (bubble = ZCTA population, colour = smoke days),
which only checks the pipeline; the proxy figure is NOT for the manuscript."""
import sys, os, numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "local_exposure_prep", "raw")
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 9})

# zip3 population-weighted centroid from ZCTA gazetteer
gaz = pd.read_csv(os.path.join(RAW, "2019_Gaz_zcta_national.txt"), sep="\t", dtype={"GEOID": str})
gaz.columns = [c.strip() for c in gaz.columns]
pop = pd.read_csv(os.path.join(RAW, "zcta_population.csv"), dtype={"zcta": str})
pop["zcta"] = pop["zcta"].str.zfill(5)
gaz = gaz.merge(pop, left_on="GEOID", right_on="zcta", how="left").fillna({"pop": 0})
gaz["zip3"] = gaz["GEOID"].str[:3]; gaz["w"] = gaz["pop"].clip(lower=1)
cen = gaz.groupby("zip3").apply(lambda d: pd.Series({"lon": np.average(d["INTPTLONG"], weights=d["w"]),
                                                     "lat": np.average(d["INTPTLAT"], weights=d["w"]),
                                                     "pop": d["pop"].sum()})).reset_index()

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results_import", "rev6_zip3_map.csv")
proxy = not os.path.exists(src)
if proxy:
    import pyarrow.parquet as pq
    s = pq.read_table(os.path.join(ROOT, "local_exposure_prep", "out", "smoke_zip3_daily.parquet")).to_pandas()
    s = s[(s["date"] >= "2016-01-01")]
    d = s.groupby("zip3").agg(colour=("smokePM_popwt", lambda x: (x > 20).sum())).reset_index()
    d = d.merge(cen, on="zip3"); d["size"] = d["pop"]
    size_label, colour_label, size_ticks = "ZCTA population", "Smoke days > 20 µg/m³, 2016-2023", [100000, 500000, 1000000]
    out = os.path.join(ROOT, "figures", "Fig_map_zip3_PROXY")
else:
    d = pd.read_csv(src, dtype={"zip3": str}).dropna(subset=["n_persons"])
    d = d.merge(cen, on="zip3"); d["size"] = d["n_persons"]; d["colour"] = d["tst_mean"]
    size_label, colour_label, size_ticks = "Participants", "Mean total sleep time, min", [30, 100, 300, 1000]
    out = os.path.join(ROOT, "figures", "Fig_map_zip3")

st = gpd.read_file(os.path.join(RAW, "us_states", "cb_2023_us_state_20m.shp"))
st = st[~st["STUSPS"].isin(["AK", "HI", "PR", "VI", "GU", "MP", "AS"])].to_crs("EPSG:5070")
pts = gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d["lon"], d["lat"]), crs="EPSG:4326").to_crs("EPSG:5070")
pts = pts[pts.within(st.union_all().buffer(20000))]

def scale(n): return 6 + 60 * np.sqrt(n / max(size_ticks))
fig, ax = plt.subplots(figsize=(7.2, 4.6))
st.boundary.plot(ax=ax, linewidth=0.4, color="#666666")
pts = pts.sort_values("size", ascending=False)
vmin, vmax = np.nanpercentile(pts["colour"], [5, 95])
sc = ax.scatter(pts.geometry.x, pts.geometry.y, s=scale(pts["size"]), c=pts["colour"], cmap="viridis",
                vmin=vmin, vmax=vmax, edgecolor="black", linewidth=0.3, alpha=0.9)
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01); cb.set_label(colour_label); cb.outline.set_linewidth(0.5)
handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="black",
                  markersize=np.sqrt(scale(n)), label=f"{n:,}") for n in size_ticks]
ax.legend(handles=handles, title=size_label, loc="lower left", frameon=False, labelspacing=1.4, borderpad=0.6)
ax.set_axis_off(); fig.tight_layout()
for ext in ("png", "svg", "pdf"): fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print("saved", out, "| points:", len(pts), "| proxy" if proxy else "| real")
