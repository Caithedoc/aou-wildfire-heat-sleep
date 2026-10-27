"""Figure in the style of Lu 2025 Nat Health Fig 2: national ZIP3 dots + eight metro-area ZCTA choropleths.
Quantity mapped: number of hot-and-smoky days per ZCTA, 2016-2023 (smoke PM2.5 > 20 µg/m³ and daily mean
temperature >= ZCTA x calendar-month 90th percentile). All inputs are public (local_exposure_prep/out), no AoU data.
Run: /usr/bin/python3 figures/make_fig_compound_zcta_map.py"""
import os, numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEP = os.path.join(ROOT, "local_exposure_prep"); OUT = os.path.join(LEP, "out"); RAW = os.path.join(LEP, "raw")
plt.rcParams.update({"font.family": "Helvetica", "font.size": 8})
CMAP = "RdYlGn_r"  # Lu 2025 Fig 2 palette, reversed so more hazard days = red
DAYS = pd.date_range("2016-01-01", "2023-12-31", freq="D")
CACHE = os.path.join(OUT, "zcta_compound_days_2016_2023.csv")

if not os.path.exists(CACHE):
    idx = pd.read_csv(os.path.join(OUT, "zcta_matrix_index.csv"), dtype={"zcta": str, "zip3": str})
    S = np.load(os.path.join(OUT, "zcta_smoke_matrix.npy"), mmap_mode="r")
    T = np.load(os.path.join(OUT, "zcta_tmean_matrix.npy"), mmap_mode="r")
    assert S.shape == T.shape == (len(idx), len(DAYS))
    months = DAYS.month.values
    hot = np.zeros(S.shape, dtype=bool)
    for m in range(1, 13):
        cols = np.where(months == m)[0]
        Tm = np.asarray(T[:, cols], dtype=np.float32)
        q = np.nanpercentile(Tm, 90, axis=1, keepdims=True)
        hot[:, cols] = Tm >= q
    smoke = np.asarray(S) > 20
    idx["smoke_days"] = smoke.sum(1); idx["hot_days"] = hot.sum(1); idx["compound_days"] = (smoke & hot).sum(1)
    idx["compound_days_2023"] = (smoke & hot)[:, DAYS.year.values == 2023].sum(1)
    idx.to_csv(CACHE, index=False); print("cached", CACHE)
d = pd.read_csv(CACHE, dtype={"zcta": str, "zip3": str})
print(d[["smoke_days", "hot_days", "compound_days"]].describe().round(1))

# national zip3 dots (population-weighted mean compound days) at population-weighted centroids
gaz = pd.read_csv(os.path.join(RAW, "2019_Gaz_zcta_national.txt"), sep="\t", dtype={"GEOID": str})
gaz.columns = [c.strip() for c in gaz.columns]
d = d.merge(gaz[["GEOID", "INTPTLAT", "INTPTLONG"]], left_on="zcta", right_on="GEOID", how="left")
d["w"] = d["pop"].clip(lower=1)
z3 = d.groupby("zip3").apply(lambda g: pd.Series({
    "lon": np.average(g["INTPTLONG"], weights=g["w"]), "lat": np.average(g["INTPTLAT"], weights=g["w"]),
    "compound": np.average(g["compound_days"], weights=g["w"]), "pop": g["pop"].sum()})).reset_index()

st = gpd.read_file(os.path.join(RAW, "us_states", "cb_2023_us_state_20m.shp"))
st = st[~st["STUSPS"].isin(["AK", "HI", "PR", "VI", "GU", "MP", "AS"])].to_crs("EPSG:5070")
zc = gpd.read_file(os.path.join(RAW, "us_zcta", "cb_2020_us_zcta520_500k.shp"))[["ZCTA5CE20", "geometry"]]
zc = zc.merge(d[["zcta", "compound_days"]], left_on="ZCTA5CE20", right_on="zcta", how="inner").to_crs("EPSG:5070")
cbsa = gpd.read_file(os.path.join(RAW, "us_cbsa", "cb_2023_us_cbsa_500k.shp")).to_crs("EPSG:5070")
METROS = [("Seattle", "Seattle-Tacoma"), ("Portland", "Portland-Vancouver"), ("San Francisco", "San Francisco-Oakland"),
          ("Los Angeles", "Los Angeles-Long Beach"), ("Denver", "Denver-Aurora"), ("Chicago", "Chicago-Naperville"),
          ("New York", "New York-Newark"), ("Washington, DC", "Washington-Arlington")]
vmax = float(np.nanpercentile(zc["compound_days"], 98)); cmap = CMAP

fig = plt.figure(figsize=(6.5, 7.0))
gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1, 2.1, 1], hspace=0.25, wspace=0.08)
def metro_panel(ax, label, prefix):
    poly = cbsa[cbsa["NAME"].str.startswith(prefix)]
    assert len(poly) == 1, (prefix, poly["NAME"].tolist())
    geom = poly.geometry.iloc[0]
    sub = gpd.clip(zc[zc.intersects(geom)], geom).explode(index_parts=False)
    sub = sub[sub.geometry.geom_type.isin(["Polygon", "MultiPolygon"])]
    gpd.GeoSeries([geom], crs=cbsa.crs).plot(ax=ax, color="#e6e6e6", linewidth=0)
    sub.plot(ax=ax, column="compound_days", cmap=cmap, vmin=0, vmax=vmax, linewidth=0.05, edgecolor="white")
    poly.boundary.plot(ax=ax, color="black", linewidth=0.5)
    ax.set_title(label, fontsize=8, pad=2); ax.set_axis_off()
    return poly.geometry.centroid.iloc[0]
cents = {}
for k, (label, prefix) in enumerate(METROS[:4]): cents[label] = metro_panel(fig.add_subplot(gs[0, k]), label, prefix)
for k, (label, prefix) in enumerate(METROS[4:]): cents[label] = metro_panel(fig.add_subplot(gs[2, k]), label, prefix)
axm = fig.add_subplot(gs[1, :])
st.boundary.plot(ax=axm, linewidth=0.35, color="#777777")
pts = gpd.GeoDataFrame(z3, geometry=gpd.points_from_xy(z3["lon"], z3["lat"]), crs="EPSG:4326").to_crs("EPSG:5070")
pts = pts[pts.within(st.union_all().buffer(20000))]
sc = axm.scatter(pts.geometry.x, pts.geometry.y, s=9, c=pts["compound"], cmap=cmap, vmin=0, vmax=vmax,
                 edgecolor="black", linewidth=0.15)
for label, c in cents.items():
    axm.scatter([c.x], [c.y], s=40, facecolor="none", edgecolor="black", linewidth=0.9)
    axm.annotate(label, (c.x, c.y), xytext=(4, 4), textcoords="offset points", fontsize=6.5)
axm.set_axis_off()
cax = fig.add_axes([0.905, 0.36, 0.014, 0.28]); cb = fig.colorbar(sc, cax=cax); cb.set_label("Hot-and-smoky days, 2016-2023", fontsize=7.5)
cb.set_ticks([0, 5, 10, 15, 20])
cb.ax.tick_params(labelsize=7); cb.outline.set_linewidth(0.5)
out = os.path.join(ROOT, "figures", "Fig_compound_days_zcta_map")
for ext in ("png", "svg", "pdf"): fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
print("saved", out)
