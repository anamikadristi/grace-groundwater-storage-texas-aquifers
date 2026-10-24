"""Build aquifer regions and area weights on the GRACE (0.25 deg) and CPC (0.5 deg) grids.

For each of the four major aquifers (Texas portion only, outcrop plus
subsurface, as mapped by TWDB), the script computes the fraction of every grid
cell that falls inside the aquifer, in an equal-area projection. Weights are
fraction x cell area, so an aquifer mean is an area-weighted mean.
"""
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
GIS = ROOT / "outputs" / "gis"
PROC.mkdir(parents=True, exist_ok=True)
GIS.mkdir(parents=True, exist_ok=True)

EA = "EPSG:3083"  # NAD83 Texas Centric Albers Equal Area
AQ = {
    "OGALLALA": "Ogallala (High Plains)",
    "GULF_COAST": "Gulf Coast",
    "EDWARDS-TRINITY": "Edwards-Trinity (Plateau)",
    "CARRIZO": "Carrizo-Wilcox",
}


def load_regions():
    g = gpd.read_file(RAW / "aquifers" / "NEW_major_aquifers_dd.shp")
    g = g[g["AQ_NAME"].isin(AQ) & g["AQUIFER"].isin([1, 2])]
    g = g.dissolve("AQ_NAME").reset_index()
    tx = gpd.read_file(RAW / "census" / "cb_2023_us_state_500k.shp")
    tx = tx[tx["STUSPS"] == "TX"].to_crs(g.crs)
    g = gpd.clip(g, tx)
    g["aquifer"] = g["AQ_NAME"].map(AQ)
    g = g.to_crs(EA)
    g["geometry"] = g.buffer(0)
    g["area_km2"] = g.area / 1e6
    return g[["aquifer", "area_km2", "geometry"]], tx.to_crs(EA)


def cell_weights(regions, res, lat_rng=(25.0, 37.0), lon_rng=(-107.0, -93.0)):
    lats = np.arange(lat_rng[0] + res / 2, lat_rng[1], res)
    lons = np.arange(lon_rng[0] + res / 2, lon_rng[1], res)
    cells = [box(lo - res / 2, la - res / 2, lo + res / 2, la + res / 2)
             for la in lats for lo in lons]
    grid = gpd.GeoDataFrame(
        {"lat": np.repeat(lats, len(lons)), "lon": np.tile(lons, len(lats))},
        geometry=cells, crs="EPSG:4326").to_crs(EA)
    grid["cell_km2"] = grid.area / 1e6
    rows = []
    for _, r in regions.iterrows():
        sub = grid[grid.intersects(r.geometry)].copy()
        sub["in_km2"] = sub.intersection(r.geometry).area / 1e6
        sub["aquifer"] = r.aquifer
        rows.append(sub[["aquifer", "lat", "lon", "cell_km2", "in_km2"]])
    w = pd.concat(rows, ignore_index=True)
    w["frac"] = w["in_km2"] / w["cell_km2"]
    return w[w["in_km2"] > 0]


def main():
    regions, tx = load_regions()
    regions.to_crs("EPSG:4326").to_file(GIS / "aquifer_regions.gpkg", layer="aquifers", driver="GPKG")
    tx.to_crs("EPSG:4326").to_file(GIS / "aquifer_regions.gpkg", layer="texas", driver="GPKG")
    print(regions[["aquifer", "area_km2"]].round(0))
    for res, name in [(0.25, "grace"), (0.5, "cpc")]:
        w = cell_weights(regions, res)
        w.to_csv(PROC / f"weights_{name}_grid.csv", index=False)
        print(name, w.groupby("aquifer")["in_km2"].sum().round(0).to_dict())
    regions.drop(columns="geometry").to_csv(PROC / "aquifer_areas.csv", index=False)


if __name__ == "__main__":
    main()
