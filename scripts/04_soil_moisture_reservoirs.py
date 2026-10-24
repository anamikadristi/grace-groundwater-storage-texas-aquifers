"""Soil moisture and reservoir storage anomalies for each aquifer region, and the
groundwater storage anomaly (GWSA) by subtraction:

    GWSA = TWSA - SMA - RSA        (snow is neglected, see README)

Soil moisture: NOAA CPC monthly soil moisture v2 (leaky bucket model, one
1.6 m column, 0.5 deg). Reservoirs: TWDB Water Data for Texas daily storage
for the 122 major reservoirs, monthly mean, summed over reservoirs located in
each aquifer region and spread over the region area.

All anomalies use the same 2004.0 to 2009.999 baseline as the CSR mascons.
"""
import json
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
TAB = ROOT / "outputs" / "tables"
GIS = ROOT / "outputs" / "gis"
EA = "EPSG:3083"
ACFT_M3 = 1233.48184
BASE = ("2004-01", "2009-12")


def soil_moisture():
    sm = xr.open_dataset(RAW / "cpc" / "cpc_soilw_texas.nc")["soilw"].load()
    w = pd.read_csv(PROC / "weights_cpc_grid.csv")
    w["lon360"] = w["lon"] % 360
    months = pd.PeriodIndex(pd.to_datetime(sm.time.values), freq="M")
    out = []
    for aq, g in w.groupby("aquifer"):
        v = sm.sel(lat=xr.DataArray(g["lat"].values, dims="c"),
                   lon=xr.DataArray(g["lon360"].values, dims="c"), method="nearest").values
        ok = np.isfinite(v).all(axis=0)
        wt = g["in_km2"].values * ok
        ts = np.nansum(v * wt, axis=1) / wt.sum()
        out.append(pd.DataFrame({"month": months.astype(str), "aquifer": aq, "sm_mm": ts}))
    df = pd.concat(out)
    df = df[df.month >= "2002-01"]
    base = df[(df.month >= BASE[0]) & (df.month <= BASE[1])].groupby("aquifer")["sm_mm"].mean()
    df["sma_mm"] = df["sm_mm"] - df["aquifer"].map(base)
    return df


def load_reservoir_points():
    d = RAW / "reservoirs"
    gj = json.loads((d / "recent-conditions.geojson").read_text())
    pts = gpd.GeoDataFrame(
        [{"short_name": f["properties"]["short_name"], "full_name": f["properties"]["full_name"],
          "conservation_capacity_acft": f["properties"]["conservation_capacity"]}
         for f in gj["features"]],
        geometry=gpd.points_from_xy([f["geometry"]["coordinates"][0] for f in gj["features"]],
                                    [f["geometry"]["coordinates"][1] for f in gj["features"]]),
        crs="EPSG:4326")
    # slug used by the CSV file names
    pts["slug"] = (pts.short_name.str.lower().str.replace("'", "", regex=False)
                   .str.replace(r"[^a-z0-9]+", "-", regex=True).str.strip("-"))
    series, meta = {}, []
    for _, r in pts.iterrows():
        f = d / "individual" / f"{r.slug}.csv"
        if not f.exists():
            meta.append({"reservoir": r.full_name, "slug": r.slug, "file": False})
            continue
        x = pd.read_csv(f, comment="#", usecols=["date", "reservoir_storage"], parse_dates=["date"]).dropna()
        m = x.set_index("date")["reservoir_storage"].resample("MS").mean()
        m.index = m.index.to_period("M")
        series[r.slug] = m
        meta.append({"reservoir": r.full_name, "slug": r.slug, "file": True,
                     "first_date": x.date.min().date(), "last_date": x.date.max().date(),
                     "n_days": len(x), "lon": round(r.geometry.x, 5), "lat": round(r.geometry.y, 5),
                     "usable": x.date.min() <= pd.Timestamp("2004-01-31")})
    return pts.to_crs(EA), series, pd.DataFrame(meta)


def reservoir_anomaly(pts, series, regions, label):
    """Sum monthly storage of usable reservoirs located in each region.

    A reservoir is usable only if its TWDB record starts by January 2004, so it
    has a full 2004 to 2009 baseline. Reservoirs whose record starts later are
    left out rather than assumed empty."""
    j = gpd.sjoin(pts, regions[["aquifer", "geometry"]], how="inner", predicate="within")
    idx = pd.period_range("2002-01", "2026-12", freq="M")
    areas = pd.read_csv(PROC / "aquifer_areas.csv").set_index("aquifer")["area_km2"]
    out, memb = [], []
    for aq in areas.index:
        slugs = [s for s in j.loc[j.aquifer == aq, "slug"] if s in series
                 and series[s].index.min() <= pd.Period("2004-01", "M")]
        memb += [{"slug": s, "region": aq, "set": label} for s in slugs]
        tot = pd.Series(0.0, index=idx)
        for s in slugs:
            m = series[s].reindex(idx).interpolate(limit_area="inside")
            # GRACE starts 2002-04; fill any month before a record start of
            # 2002-2003 with the first value (affects at most a few months)
            m = m.bfill()
            tot = tot + m
        last = min(series[s].index.max() for s in slugs) if slugs else idx[-1]
        tot[tot.index > last] = np.nan
        mm = tot * ACFT_M3 / (areas[aq] * 1e6) * 1000.0
        out.append(pd.DataFrame({"month": idx.astype(str), "aquifer": aq,
                                 f"res_storage_acft_{label}": tot.values, f"res_mm_{label}": mm.values,
                                 f"n_res_{label}": len(slugs)}))
    df = pd.concat(out)
    base = df[(df.month >= BASE[0]) & (df.month <= BASE[1])].groupby("aquifer")[f"res_mm_{label}"].mean()
    df[f"rsa_mm_{label}"] = df[f"res_mm_{label}"] - df["aquifer"].map(base)
    return df, pd.DataFrame(memb)


def main():
    sm = soil_moisture()
    pts, series, meta = load_reservoir_points()
    regions = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="aquifers").to_crs(EA)
    buf = regions.copy()
    buf["geometry"] = buf.buffer(100_000)
    res_in, m_in = reservoir_anomaly(pts, series, regions, "in")
    res_buf, m_buf = reservoir_anomaly(pts, series, buf, "buf100")
    memb = pd.concat([m_in, m_buf])
    memb = memb.pivot_table(index="slug", columns="set", values="region",
                            aggfunc=lambda s: "; ".join(sorted(s))).reset_index()
    meta = meta.merge(memb, on="slug", how="left").rename(
        columns={"in": "regions_inside", "buf100": "regions_within_100km"})
    meta.to_csv(TAB / "reservoirs_used.csv", index=False)
    print("reservoirs with a record starting after Jan 2004 (left out):",
          meta.loc[meta.usable == False, "reservoir"].tolist())
    print("usable reservoirs inside each aquifer:", res_in.groupby("aquifer").n_res_in.first().to_dict())
    print("usable reservoirs within 100 km:", res_buf.groupby("aquifer").n_res_buf100.first().to_dict())

    tw = pd.read_csv(PROC / "grace_twsa_aquifer_monthly.csv")
    df = tw.merge(sm[["month", "aquifer", "sm_mm", "sma_mm"]], on=["month", "aquifer"], how="left")
    df = df.merge(res_in.drop(columns="n_res_in"), on=["month", "aquifer"], how="left")
    df = df.merge(res_buf.drop(columns="n_res_buf100"), on=["month", "aquifer"], how="left")
    df["gwsa_mm"] = df["twsa_mm"] - df["sma_mm"] - df["rsa_mm_in"]
    df["gwsa_mm_resbuf100"] = df["twsa_mm"] - df["sma_mm"] - df["rsa_mm_buf100"]
    df.to_csv(PROC / "storage_components_aquifer_monthly.csv", index=False)
    print(df.groupby("aquifer")[["twsa_mm", "sma_mm", "rsa_mm_in", "rsa_mm_buf100", "gwsa_mm"]].std().round(1))
    print("months with GWSA:", df.gwsa_mm.notna().groupby(df.aquifer).sum().to_dict())


if __name__ == "__main__":
    main()
