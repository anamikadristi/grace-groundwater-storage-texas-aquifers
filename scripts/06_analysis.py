"""Trends, STL decomposition, drought response, well correlation and a
resolution and leakage check.

Outputs (outputs/tables):
  table1_trends.csv               linear trends by component, period, aquifer
  table2_stl_summary.csv          seasonal amplitude and STL trend change
  table3_drought_response.csv     2011 and 2022 drought drawdown and recovery
  table4_well_correlation.csv     GRACE vs TWDB well-level anomalies
  table5_leakage_sensitivity.csv  aquifer masks shrunk and grown
  table6_resolution.csv           aquifer size vs GRACE resolution, overlaps
  data/processed/stl_components_monthly.csv
"""
import importlib.util
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr
from scipy import stats
from statsmodels.tsa.seasonal import STL

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
TAB = ROOT / "outputs" / "tables"
GIS = ROOT / "outputs" / "gis"
EA = "EPSG:3083"
PERIODS = {
    "full 2002-04 to 2026-08": ("2002-04", "2026-08"),
    "GRACE 2002-04 to 2017-06": ("2002-04", "2017-06"),
    "GRACE-FO 2018-06 to 2026-08": ("2018-06", "2026-08"),
}
COMPS = {"twsa_mm": "TWSA", "sma_mm": "Soil moisture", "rsa_mm_in": "Reservoirs",
         "gwsa_mm": "GWSA", "gwsa_mm_resbuf100": "GWSA (reservoirs within 100 km)"}


def tyear(month_str):
    p = pd.PeriodIndex(month_str, freq="M")
    return np.asarray(p.year + (p.month - 0.5) / 12.0, dtype=float)


def trend(t, y, maxlags=12):
    ok = np.isfinite(y)
    t, y = t[ok], y[ok]
    X = np.column_stack([np.ones_like(t), t - 2014.0, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t),
                         np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)])
    f = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": maxlags})
    lo, hi = f.conf_int()[1]
    return dict(n_months=int(ok.sum()), trend_mm_yr=f.params[1], se_mm_yr=f.bse[1],
                ci95_lo=lo, ci95_hi=hi, p_value=f.pvalues[1])


def regular(df, col):
    """Regular monthly series, linear interpolation across missing months."""
    s = df.set_index(pd.PeriodIndex(df.month, freq="M"))[col]
    idx = pd.period_range(s.index.min(), s.index.max(), freq="M")
    s = s.reindex(idx)
    return s, s.interpolate(limit_area="inside")


def main():
    df = pd.read_csv(PROC / "storage_components_aquifer_monthly.csv")
    areas = pd.read_csv(PROC / "aquifer_areas.csv").set_index("aquifer")["area_km2"]
    aquifers = list(areas.index)

    # ---------------- 1. trends ----------------
    rows = []
    for aq in aquifers:
        d = df[df.aquifer == aq]
        for pname, (a, b) in PERIODS.items():
            dd = d[(d.month >= a) & (d.month <= b)]
            t = tyear(dd.month)
            for c, cname in COMPS.items():
                r = trend(t, dd[c].values)
                r.update(aquifer=aq, period=pname, component=cname,
                         trend_km3_yr=r["trend_mm_yr"] * areas[aq] * 1e-6)
                rows.append(r)
    t1 = pd.DataFrame(rows)[["aquifer", "component", "period", "n_months", "trend_mm_yr", "se_mm_yr",
                             "ci95_lo", "ci95_hi", "p_value", "trend_km3_yr"]]
    t1.to_csv(TAB / "table1_trends.csv", index=False, float_format="%.4g")

    # ---------------- 2. STL ----------------
    stl_rows, stl_out = [], []
    for aq in aquifers:
        d = df[df.aquifer == aq]
        for c in ["twsa_mm", "gwsa_mm"]:
            raw, filled = regular(d, c)
            res = STL(filled.values, period=12, robust=True, seasonal=13).fit()
            o = pd.DataFrame({"month": filled.index.astype(str), "aquifer": aq, "component": COMPS[c],
                              "observed": raw.values, "filled": filled.values, "trend": res.trend,
                              "seasonal": res.seasonal, "resid": res.resid,
                              "interpolated": raw.isna().values})
            stl_out.append(o)
            clim = o.assign(m=filled.index.month).groupby("m").seasonal.mean()
            tr = pd.Series(res.trend, index=filled.index)
            obs_var = np.nanvar(o.filled)
            stl_rows.append({
                "aquifer": aq, "component": COMPS[c],
                "seasonal_amplitude_mm": clim.max() - clim.min(),
                "seasonal_peak_month": int(clim.idxmax()), "seasonal_trough_month": int(clim.idxmin()),
                "stl_trend_2003_mean_mm": tr[tr.index.year == 2003].mean(),
                "stl_trend_2025_mean_mm": tr[tr.index.year == 2025].mean(),
                "share_var_trend": np.var(res.trend) / obs_var,
                "share_var_seasonal": np.var(res.seasonal) / obs_var,
                "share_var_resid": np.var(res.resid) / obs_var,
                "months_interpolated": int(raw.isna().sum())})
    stl_df = pd.concat(stl_out)
    stl_df.to_csv(PROC / "stl_components_monthly.csv", index=False, float_format="%.3f")
    pd.DataFrame(stl_rows).to_csv(TAB / "table2_stl_summary.csv", index=False, float_format="%.3f")

    # ---------------- 3. droughts ----------------
    # deseasonalised series (observed months only, minus STL seasonal) so a
    # summer low is not mistaken for drought loss
    dro = []
    events = {"2011": ("2010-01", "2010-12", "2011-01", "2012-12", "2013-01", "2013-12"),
              "2022": ("2021-01", "2021-12", "2022-01", "2023-12", "2024-01", "2024-12")}
    for aq in aquifers:
        for comp in ["TWSA", "GWSA"]:
            o = stl_df[(stl_df.aquifer == aq) & (stl_df.component == comp) & ~stl_df.interpolated].copy()
            o["deseas"] = o.observed - o.seasonal
            for ev, (r0, r1, w0, w1, a0, a1) in events.items():
                ref = o[(o.month >= r0) & (o.month <= r1)].deseas.mean()
                win = o[(o.month >= w0) & (o.month <= w1)]
                imin = win.deseas.idxmin()
                after = o[(o.month >= a0) & (o.month <= a1)].deseas.mean()
                drop = win.loc[imin, "deseas"] - ref
                dro.append({"aquifer": aq, "component": comp, "drought": ev,
                            "reference_period": f"{r0} to {r1}", "reference_mm": ref,
                            "min_month": win.loc[imin, "month"], "min_mm": win.loc[imin, "deseas"],
                            "drawdown_mm": drop, "drawdown_km3": drop * areas[aq] * 1e-6,
                            "following_year_mean_mm": after,
                            # recovery share only when the drawdown is larger than 20 mm
                            "recovered_fraction": (after - win.loc[imin, "deseas"]) / -drop if drop < -20 else np.nan,
                            "n_months_reference": int(o[(o.month >= r0) & (o.month <= r1)].shape[0])})
    pd.DataFrame(dro).to_csv(TAB / "table3_drought_response.csv", index=False, float_format="%.3f")

    # ---------------- 4. wells ----------------
    wa = pd.read_csv(PROC / "wells_annual_composite.csv")
    wm = pd.read_csv(PROC / "wells_monthly_composite.csv")
    df["year"] = df.month.str[:4].astype(int)
    ann = df.groupby(["aquifer", "year"]).agg(n=("gwsa_mm", "count"), gwsa_mm=("gwsa_mm", "mean"),
                                              twsa_mm=("twsa_mm", "mean")).reset_index()
    ann = ann[ann.n >= 6]
    rows = []
    for aq in aquifers:
        m = ann[ann.aquifer == aq].merge(wa[wa.aquifer == aq], on=["aquifer", "year"]).sort_values("year")
        for comp in ["gwsa_mm", "twsa_mm"]:
            x, y = m.well_anom_m.values, m[comp].values
            r, p = stats.pearsonr(x, y)
            rho, prho = stats.spearmanr(x, y)
            dx, dy = np.diff(x), np.diff(y)
            consec = np.diff(m.year.values) == 1
            rd, pd_ = stats.pearsonr(dx[consec], dy[consec])
            slope = stats.linregress(x, y).slope
            rows.append({"aquifer": aq, "grace_component": comp.split("_")[0].upper(), "scale": "annual",
                         "n": len(m), "years": f"{m.year.min()}-{m.year.max()}", "pearson_r": r, "p": p,
                         "spearman_rho": rho, "p_spearman": prho, "r_first_differences": rd, "p_first_diff": pd_,
                         "slope_mm_per_m": slope, "apparent_storage_coefficient": slope / 1000.0,
                         "median_wells_per_year": float(m.n_wells.median())})
        mm = df[df.aquifer == aq].merge(wm[(wm.aquifer == aq) & (wm.n_wells >= 5)], on=["aquifer", "month"])
        for comp in ["gwsa_mm", "twsa_mm"]:
            x, y = mm.well_anom_m.values, mm[comp].values
            if len(mm) < 24:
                continue
            r, p = stats.pearsonr(x, y)
            rho, prho = stats.spearmanr(x, y)
            rows.append({"aquifer": aq, "grace_component": comp.split("_")[0].upper(), "scale": "monthly",
                         "n": len(mm), "years": f"{mm.month.min()} to {mm.month.max()}", "pearson_r": r, "p": p,
                         "spearman_rho": rho, "p_spearman": prho,
                         "slope_mm_per_m": stats.linregress(x, y).slope,
                         "apparent_storage_coefficient": stats.linregress(x, y).slope / 1000.0,
                         "median_wells_per_year": float(mm.n_wells.median())})
    pd.DataFrame(rows).to_csv(TAB / "table4_well_correlation.csv", index=False, float_format="%.3f")
    ann.merge(wa, on=["aquifer", "year"], how="left").to_csv(PROC / "annual_grace_vs_wells.csv", index=False,
                                                              float_format="%.4f")

    # ---------------- 5. leakage sensitivity ----------------
    spec = importlib.util.spec_from_file_location("regions", Path(__file__).with_name("02_aquifer_regions.py"))
    regmod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(regmod)
    reg = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="aquifers").to_crs(EA)
    ds = xr.open_dataset(RAW / "grace" / "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc", decode_times=False)
    sub = ds["lwe_thickness"].sel(lat=slice(24, 39), lon=slice(250, 269)).load() * 10.0
    land = xr.open_dataset(RAW / "grace" / "CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc")["LO_val"]
    land = land.sel(lat=slice(24, 39), lon=slice(250, 269)).load()
    months = pd.read_csv(PROC / "grace_twsa_aquifer_monthly.csv").drop_duplicates("month").month.values
    t = tyear(months)
    base = df.set_index(["aquifer", "month"])
    rows = []
    for dist in [-50, 0, 100, 200]:
        r2 = reg.copy()
        r2["geometry"] = r2.buffer(dist * 1000) if dist else r2.geometry
        r2 = r2[~r2.geometry.is_empty]
        w = regmod.cell_weights(r2, 0.25, lat_rng=(24.0, 39.0), lon_rng=(-110.0, -91.0))
        for aq, g in w.groupby("aquifer"):
            sel = dict(lat=xr.DataArray(g.lat.values, dims="c"), lon=xr.DataArray(g.lon.values % 360, dims="c"))
            v = sub.sel(**sel, method="nearest").values
            lm = land.sel(**sel, method="nearest").values == 1
            wt = g.in_km2.values * lm
            tw = (v * wt).sum(axis=1) / wt.sum()
            gws = tw - base.loc[aq].loc[months, "sma_mm"].values - base.loc[aq].loc[months, "rsa_mm_in"].values
            r_tw = trend(t, tw)
            r_gw = trend(t, gws)
            rows.append({"aquifer": aq, "mask_buffer_km": dist, "mask_area_km2": float(g.in_km2.sum()),
                         "twsa_trend_mm_yr": r_tw["trend_mm_yr"], "twsa_se": r_tw["se_mm_yr"],
                         "gwsa_trend_mm_yr": r_gw["trend_mm_yr"], "gwsa_se": r_gw["se_mm_yr"],
                         "r_with_base_twsa": np.corrcoef(tw, base.loc[aq].loc[months, "twsa_mm"].values)[0, 1]})
    t5 = pd.DataFrame(rows).sort_values(["aquifer", "mask_buffer_km"])
    t5.to_csv(TAB / "table5_leakage_sensitivity.csv", index=False, float_format="%.3f")

    # ---------------- 6. resolution and overlap ----------------
    rows = []
    for _, r in reg.iterrows():
        others = reg[reg.aquifer != r.aquifer]
        ov = {f"overlap_with_{o.aquifer}_km2": r.geometry.intersection(o.geometry).area / 1e6
              for _, o in others.iterrows()}
        rows.append({"aquifer": r.aquifer, "area_km2": r.geometry.area / 1e6,
                     "area_over_300km_footprint": (r.geometry.area / 1e6) / (300 ** 2),
                     "max_inscribed_radius_km": _inscribed_radius_km(r.geometry), **ov})
    t6 = pd.DataFrame(rows)
    corr = df.pivot(index="month", columns="aquifer", values="twsa_mm").corr()
    t6 = t6.merge(corr.add_prefix("twsa_r_with_").reset_index().rename(columns={"index": "aquifer"}),
                  on="aquifer")
    t6.to_csv(TAB / "table6_resolution.csv", index=False, float_format="%.3f")
    print(t1[(t1.component.isin(["TWSA", "GWSA"]))].round(3).to_string())


def _inscribed_radius_km(geom):
    """Largest distance from any interior point to the boundary (approximate)."""
    from shapely.ops import polylabel
    g = max(geom.geoms, key=lambda x: x.area) if geom.geom_type == "MultiPolygon" else geom
    p = polylabel(g, tolerance=1000)
    return g.boundary.distance(p) / 1000


if __name__ == "__main__":
    main()
