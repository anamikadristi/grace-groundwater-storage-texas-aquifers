"""In-situ groundwater level anomalies by aquifer from the TWDB Groundwater
Database (GWDB), for comparison with GRACE GWSA.

Screening:
  * Status "Publishable", numeric depth to water, 2002 onward
  * pumping or recently pumped measurements removed (Remarks)
  * well's GWDB aquifer is one of the four study aquifers and the well location
    falls inside that aquifer region (Texas part)

Two composites per aquifer:
  annual  : wells with an annual median in at least 15 years of 2002 to 2025
            and at least 3 years in the 2004 to 2009 baseline. Anomaly = annual
            median water level minus the well's 2004 to 2009 mean. Aquifer
            value = median across wells.
  monthly : wells with at least 120 monthly values from 2002 onward (mostly
            recorder and transducer wells). Anomaly = monthly median minus the
            well's mean over 2004 to 2009 if available, else over its record,
            then re-centred on the composite's 2004 to 2009 mean.
Water level = minus depth to water, converted from feet to metres.
"""
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
TAB = ROOT / "outputs" / "tables"
GIS = ROOT / "outputs" / "gis"
FT = 0.3048
AQMAP = {
    "Ogallala": "Ogallala (High Plains)",
    "Gulf Coast": "Gulf Coast",
    "Edwards-Trinity Plateau": "Edwards-Trinity (Plateau)",
    "Carrizo-Wilcox": "Carrizo-Wilcox",
}
PUMP = r"(?i)pumping|pumped recently|pumping nearby|wells pumping"


def main():
    wl = pd.read_csv(RAW / "gwdb" / "WaterLevelsMajor.txt", sep="|", dtype=str, encoding="latin-1",
                     usecols=["StateWellNumber", "Aquifer", "Status", "MeasurementDate",
                              "MeasurementYear", "DepthFromLSD", "Remarks", "MethodOfMeasurement"])
    n0 = len(wl)
    wl = wl[wl.Aquifer.isin(AQMAP)]
    wl["year"] = pd.to_numeric(wl.MeasurementYear, errors="coerce")
    wl = wl[wl.year >= 2002]
    n1 = len(wl)
    wl = wl[wl.Status == "Publishable"]
    wl = wl[~wl.Remarks.fillna("").str.contains(PUMP, regex=True)]
    wl["depth_ft"] = pd.to_numeric(wl.DepthFromLSD, errors="coerce")
    wl["date"] = pd.to_datetime(wl.MeasurementDate, errors="coerce")
    wl = wl.dropna(subset=["depth_ft", "date"])
    wl = wl[wl.date <= pd.Timestamp("2026-10-31")]
    wl["aquifer"] = wl.Aquifer.map(AQMAP)
    wl["level_m"] = -wl.depth_ft * FT
    n2 = len(wl)

    wm = pd.read_csv(RAW / "gwdb" / "WellMain.txt", sep="|", dtype=str, encoding="latin-1",
                     usecols=["StateWellNumber", "LatitudeDD", "LongitudeDD", "WellDepth", "WellUse", "County"])
    wm["lat"] = pd.to_numeric(wm.LatitudeDD, errors="coerce")
    wm["lon"] = pd.to_numeric(wm.LongitudeDD, errors="coerce")
    wells = wl[["StateWellNumber", "aquifer"]].drop_duplicates().merge(wm, on="StateWellNumber", how="left")
    wells = wells.dropna(subset=["lat", "lon"])
    g = gpd.GeoDataFrame(wells, geometry=gpd.points_from_xy(wells.lon, wells.lat), crs="EPSG:4326")
    reg = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="aquifers").to_crs("EPSG:4326")
    j = gpd.sjoin(g, reg[["aquifer", "geometry"]].rename(columns={"aquifer": "reg"}), predicate="within")
    j = j[j.aquifer == j.reg]
    keep = set(zip(j.StateWellNumber, j.aquifer))
    wl = wl[[k in keep for k in zip(wl.StateWellNumber, wl.aquifer)]]
    n3 = len(wl)
    print(f"records: all major {n0:,}; study aquifers 2002+ {n1:,}; after screening {n2:,}; inside region {n3:,}")

    # ---- annual composite ---------------------------------------------
    ann = wl.groupby(["aquifer", "StateWellNumber", "year"])["level_m"].median().reset_index()
    ann = ann[ann.year <= 2025]
    nyr = ann.groupby(["aquifer", "StateWellNumber"]).year.nunique()
    nbase = ann[ann.year.between(2004, 2009)].groupby(["aquifer", "StateWellNumber"]).year.nunique()
    ok = (nyr >= 15) & (nbase.reindex(nyr.index).fillna(0) >= 3)
    sel = ok[ok].index
    ann = ann.set_index(["aquifer", "StateWellNumber"]).loc[sel].reset_index()
    base = ann[ann.year.between(2004, 2009)].groupby(["aquifer", "StateWellNumber"]).level_m.mean()
    ann["anom_m"] = ann.level_m - pd.Series(list(zip(ann.aquifer, ann.StateWellNumber))).map(base).values
    comp = ann.groupby(["aquifer", "year"]).agg(
        well_anom_m=("anom_m", "median"), q25=("anom_m", lambda s: s.quantile(.25)),
        q75=("anom_m", lambda s: s.quantile(.75)), n_wells=("anom_m", "size")).reset_index()
    comp.to_csv(PROC / "wells_annual_composite.csv", index=False)

    # well list for maps (annual set)
    wl_sel = j.set_index(["aquifer", "StateWellNumber"]).loc[sel].reset_index()
    wl_sel = wl_sel[["StateWellNumber", "aquifer", "County", "lat", "lon", "WellDepth", "WellUse"]]
    wl_sel = wl_sel.merge(nyr.rename("n_years").reset_index(), on=["aquifer", "StateWellNumber"])
    wl_sel.to_csv(PROC / "wells_annual_set.csv", index=False)

    # ---- monthly composite --------------------------------------------
    wl["month"] = wl.date.dt.to_period("M")
    mon = wl.groupby(["aquifer", "StateWellNumber", "month"])["level_m"].median().reset_index()
    nmon = mon.groupby(["aquifer", "StateWellNumber"]).month.nunique()
    msel = nmon[nmon >= 120].index
    mon = mon.set_index(["aquifer", "StateWellNumber"]).loc[msel].reset_index()
    mon["is_base"] = (mon.month >= pd.Period("2004-01", "M")) & (mon.month <= pd.Period("2009-12", "M"))
    key = list(zip(mon.aquifer, mon.StateWellNumber))
    bmean = mon[mon.is_base].groupby(["aquifer", "StateWellNumber"]).level_m.mean()
    nb = mon[mon.is_base].groupby(["aquifer", "StateWellNumber"]).size()
    allmean = mon.groupby(["aquifer", "StateWellNumber"]).level_m.mean()
    ref = allmean.copy()
    good = nb[nb >= 24].index
    ref.loc[good] = bmean.loc[good]
    mon["anom_m"] = mon.level_m - pd.Series(key).map(ref).values
    mc = mon.groupby(["aquifer", "month"]).agg(well_anom_m=("anom_m", "median"),
                                               n_wells=("anom_m", "size")).reset_index()
    # re-centre composite on its 2004-2009 mean where available
    for aq in mc.aquifer.unique():
        m = mc.aquifer == aq
        b = mc.loc[m & (mc.month >= pd.Period("2004-01", "M")) & (mc.month <= pd.Period("2009-12", "M")), "well_anom_m"]
        if b.notna().sum() >= 24:
            mc.loc[m, "well_anom_m"] -= b.mean()
    mc["month"] = mc.month.astype(str)
    mc.to_csv(PROC / "wells_monthly_composite.csv", index=False)

    s = pd.DataFrame({
        "wells_screened_records": wl.groupby("aquifer").size(),
        "wells_with_any_record": wl.groupby("aquifer").StateWellNumber.nunique(),
        "annual_set_wells": pd.Series([a for a, _ in sel]).value_counts(),
        "monthly_set_wells": pd.Series([a for a, _ in msel]).value_counts(),
    }).fillna(0).astype(int)
    s.index.name = "aquifer"
    s.to_csv(TAB / "well_selection_summary.csv")
    print(s)


if __name__ == "__main__":
    main()
