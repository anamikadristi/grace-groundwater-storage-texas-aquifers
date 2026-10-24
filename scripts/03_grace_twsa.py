"""GRACE/GRACE-FO terrestrial water storage anomaly (TWSA) for each aquifer, and
a TWSA trend for every CSR mascon over Texas.

Input: CSR RL06.3 mascons (all corrections applied: GIA ICE6G-D, TN-14 C20/C30,
TN-13 degree 1), liquid water equivalent thickness in cm, anomaly relative to
the 2004.0 to 2009.999 mean.

Outputs:
  data/processed/grace_twsa_aquifer_monthly.csv   (mm)
  data/processed/grace_mascon_trends_texas.csv
  outputs/gis/grace_mascon_trends_texas.gpkg      (mascon polygons with trend)
  outputs/gis/grace_twsa_trend_mm_per_yr.tif      (0.25 deg GeoTIFF)
"""
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import statsmodels.api as sm
import xarray as xr
from rasterio.transform import from_origin
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
GIS = ROOT / "outputs" / "gis"
TAB = ROOT / "outputs" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

GDIR = RAW / "grace"
LAT = (25.0, 37.0)
LON = (253.0, 267.0)  # 0 to 360


def month_labels(ds):
    tb = ds["time_bounds"].values  # days since 2002-01-01
    mid = pd.Timestamp("2002-01-01") + pd.to_timedelta(tb.mean(axis=1), unit="D")
    lab = mid.to_period("M")
    # When two solutions fall in the same calendar month (a few early and
    # GRACE-FO transition months), move the later one to the next free month
    # if that month is empty, which is the usual convention.
    lab = list(lab)
    for i in range(1, len(lab)):
        if lab[i] == lab[i - 1]:
            nxt = lab[i] + 1
            if i + 1 >= len(lab) or lab[i + 1] != nxt:
                lab[i] = nxt
    return pd.PeriodIndex(lab, freq="M"), mid


def hac_trend(t_years, y, maxlags=12):
    """OLS trend with annual and semiannual harmonics; Newey-West SE."""
    ok = np.isfinite(y)
    t = t_years[ok]
    X = np.column_stack([
        np.ones(ok.sum()), t - 2014.0,
        np.sin(2 * np.pi * t), np.cos(2 * np.pi * t),
        np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)])
    fit = sm.OLS(y[ok], X).fit(cov_type="HAC", cov_kwds={"maxlags": maxlags})
    return fit.params[1], fit.bse[1], fit.pvalues[1]


def main():
    ds = xr.open_dataset(GDIR / "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc",
                         decode_times=False)
    labels, mid = month_labels(ds)
    sub = ds["lwe_thickness"].sel(lat=slice(*LAT), lon=slice(*LON)).load() * 10.0  # cm -> mm
    land = xr.open_dataset(GDIR / "CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc")["LO_val"]
    land = land.sel(lat=slice(*LAT), lon=slice(*LON)).load()
    mid_ids = xr.open_dataset(GDIR / "CSR_GRACE_GRACE-FO_RL0603_mascons_mapping_file.nc")["mascon_id"]
    mid_ids = mid_ids.sel(lat=slice(*LAT), lon=slice(*LON)).load()
    print(f"GRACE solutions: {len(labels)}, {labels[0]} to {labels[-1]}, "
          f"duplicate labels after fix: {labels.duplicated().sum()}")

    # ---- aquifer means -------------------------------------------------
    w = pd.read_csv(PROC / "weights_grace_grid.csv")
    w["lon360"] = w["lon"] % 360
    out = []
    for aq, g in w.groupby("aquifer"):
        vals = sub.sel(lat=xr.DataArray(g["lat"].values, dims="c"),
                       lon=xr.DataArray(g["lon360"].values, dims="c"), method="nearest")
        lm = land.sel(lat=xr.DataArray(g["lat"].values, dims="c"),
                      lon=xr.DataArray(g["lon360"].values, dims="c"), method="nearest").values
        wt = g["in_km2"].values * (lm == 1)
        ts = (vals.values * wt).sum(axis=1) / wt.sum()
        out.append(pd.DataFrame({"month": labels.astype(str), "mid_date": mid.date,
                                 "aquifer": aq, "twsa_mm": ts}))
        n_masc = len(np.unique(mid_ids.sel(lat=xr.DataArray(g["lat"].values, dims="c"),
                                           lon=xr.DataArray(g["lon360"].values, dims="c"),
                                           method="nearest").values))
        print(f"  {aq}: {len(g)} grid cells, {n_masc} mascons touched")
    twsa = pd.concat(out)
    twsa.to_csv(PROC / "grace_twsa_aquifer_monthly.csv", index=False)

    # ---- per-mascon trend ---------------------------------------------
    t_years = np.array([p.year + (p.month - 0.5) / 12 for p in labels])
    ids = mid_ids.values
    arr = sub.values  # time, lat, lon
    rows = []
    for mid_ in np.unique(ids[np.isfinite(ids)]):
        m = ids == mid_
        series = np.nanmean(arr[:, m], axis=1)
        spread = np.nanmax(np.nanstd(arr[:, m], axis=1))
        b, se, p = hac_trend(t_years, series)
        lat_c = float(sub.lat.values[np.where(m)[0]].mean())
        lon_c = float(sub.lon.values[np.where(m)[1]].mean()) - 360
        rows.append({"mascon_id": int(mid_), "lat_c": lat_c, "lon_c": lon_c,
                     "n_cells": int(m.sum()), "within_mascon_sd_mm": spread,
                     "trend_mm_yr": b, "se_mm_yr": se, "p_value": p})
    tr = pd.DataFrame(rows)
    tr.to_csv(PROC / "grace_mascon_trends_texas.csv", index=False)
    print(f"  mascons in Texas box: {len(tr)}, max within-mascon SD {tr.within_mascon_sd_mm.max():.3g} mm")

    # mascon polygons from 0.25 deg cells
    la, lo = np.meshgrid(sub.lat.values, sub.lon.values - 360, indexing="ij")
    cells = gpd.GeoDataFrame({"mascon_id": ids.ravel().astype(int)},
                             geometry=[box(x - .125, y - .125, x + .125, y + .125)
                                       for y, x in zip(la.ravel(), lo.ravel())],
                             crs="EPSG:4326")
    poly = cells.dissolve("mascon_id").reset_index().merge(tr, on="mascon_id")
    poly.to_file(GIS / "grace_mascon_trends_texas.gpkg", layer="mascon_twsa_trend", driver="GPKG")

    # GeoTIFF of trend on 0.25 grid, north up
    tmap = dict(zip(tr.mascon_id, tr.trend_mm_yr))
    grid = np.vectorize(lambda v: tmap.get(int(v), np.nan))(ids).astype("float32")[::-1]
    transform = from_origin(float(sub.lon.values[0] - 360 - .125), float(sub.lat.values[-1] + .125), .25, .25)
    with rasterio.open(GIS / "grace_twsa_trend_mm_per_yr.tif", "w", driver="GTiff",
                       height=grid.shape[0], width=grid.shape[1], count=1, dtype="float32",
                       crs="EPSG:4326", transform=transform, nodata=np.nan, compress="deflate") as dst:
        dst.write(grid, 1)

    # month coverage table
    cov = pd.DataFrame({"month": labels.astype(str), "mid_date": mid.date})
    allm = pd.period_range(labels[0], labels[-1], freq="M")
    missing = sorted(set(allm.astype(str)) - set(cov.month))
    pd.DataFrame({"missing_month": missing}).to_csv(TAB / "grace_missing_months.csv", index=False)
    print(f"  months present {len(cov)}, missing {len(missing)} of {len(allm)}")


if __name__ == "__main__":
    main()
