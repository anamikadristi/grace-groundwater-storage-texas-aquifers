"""Download all public, keyless input data.

Everything goes into RAW_DIR (default: data/raw, which is git-ignored).
Set the environment variable GRACE_RAW_DIR to put it somewhere else.

Sources (all keyless, no account needed):
  1. CSR RL06.3 GRACE/GRACE-FO mascons, UT Austin Center for Space Research
  2. NOAA CPC monthly soil moisture v2 (NOAA PSL), Texas subset kept
  3. TWDB major aquifer polygons
  4. TWDB Water Data for Texas reservoir storage (statewide and individual CSVs)
  5. TWDB Groundwater Database (GWDB) bulk download, water levels and well metadata
  6. U.S. Census 2023 cartographic boundary states (Texas outline for maps)

GLDAS (NASA GES DISC) and NLDAS were checked and redirect to an Earthdata
login, so they are not used. See README.
"""
import json
import os
import re
import sys
import time
import zipfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GRACE_RAW_DIR", ROOT / "data" / "raw"))
RAW.mkdir(parents=True, exist_ok=True)

CSR = "https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/"
CSR_FILES = [
    "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc",
    "CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc",
    "CSR_GRACE_GRACE-FO_RL0603_mascons_mapping_file.nc",
]
CPC_FULL = "https://downloads.psl.noaa.gov/Datasets/cpcsoil/soilw.mon.mean.v2.nc"
AQUIFERS = "https://www.twdb.texas.gov/mapping/gisdata/doc/major_aquifers.zip"
GWDB = "https://www.twdb.texas.gov/groundwater/data/GWDBDownload.zip"
STATES = "https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_state_500k.zip"
WDFT = "https://www.waterdatafortexas.org"

# Texas bounding box with margin (CPC longitudes are 0 to 360)
LAT = (24.0, 38.0)
LON360 = (252.0, 268.0)


def get(url, dest, retries=3):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  exists: {dest.name}")
        return dest
    for i in range(retries):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                tmp.rename(dest)
            print(f"  downloaded: {dest.name} ({dest.stat().st_size:,} bytes)")
            return dest
        except Exception as e:  # noqa: BLE001
            print(f"  retry {i + 1} for {url}: {e}")
            time.sleep(5)
    raise RuntimeError(f"failed: {url}")


def main():
    log = {"accessed": time.strftime("%Y-%m-%d")}

    print("1. CSR GRACE/GRACE-FO RL06.3 mascons")
    d = RAW / "grace"
    d.mkdir(exist_ok=True)
    for f in CSR_FILES:
        get(CSR + f, d / f)

    print("2. NOAA CPC soil moisture v2, Texas subset")
    d = RAW / "cpc"
    d.mkdir(exist_ok=True)
    out = d / "cpc_soilw_texas.nc"
    if not out.exists():
        import xarray as xr
        # The full global file is about 256 MB; it downloads faster than an
        # OPeNDAP subset request, so fetch it and keep only the Texas box.
        full = get(CPC_FULL, d / "soilw.mon.mean.v2.nc")
        with xr.open_dataset(full) as ds:
            da = ds["soilw"]
            lat = slice(LAT[1], LAT[0]) if da.lat[0] > da.lat[-1] else slice(*LAT)
            da.sel(lat=lat, lon=slice(*LON360)).load().to_netcdf(out)
    print(f"  saved {out.name}")

    print("3. TWDB major aquifers")
    d = RAW / "aquifers"
    d.mkdir(exist_ok=True)
    z = get(AQUIFERS, d / "major_aquifers.zip")
    with zipfile.ZipFile(z) as zf:
        zf.extractall(d)

    print("4. Census state boundaries")
    d = RAW / "census"
    d.mkdir(exist_ok=True)
    z = get(STATES, d / "cb_2023_us_state_500k.zip")
    with zipfile.ZipFile(z) as zf:
        zf.extractall(d)

    print("5. Water Data for Texas reservoirs")
    d = RAW / "reservoirs"
    d.mkdir(exist_ok=True)
    get(WDFT + "/reservoirs/statewide.csv", d / "statewide.csv")
    get(WDFT + "/reservoirs/statewide/recent-conditions.geojson", d / "recent-conditions.geojson")
    html = requests.get(WDFT + "/reservoirs/statewide", timeout=120).text
    slugs = sorted(set(re.findall(r'href="/reservoirs/individual/([a-z0-9-]+)"', html)))
    print(f"  {len(slugs)} individual reservoirs listed")
    (d / "individual").mkdir(exist_ok=True)
    # Each file is generated on request and can take a minute; fetch 4 at a time.
    from concurrent.futures import ThreadPoolExecutor

    def one(s):
        try:
            get(f"{WDFT}/reservoirs/individual/{s}.csv", d / "individual" / f"{s}.csv", retries=5)
        except RuntimeError as e:
            print(f"  skipped {s}: {e}")

    with ThreadPoolExecutor(4) as ex:
        list(ex.map(one, slugs))
    log["reservoir_slugs"] = slugs

    print("6. TWDB Groundwater Database bulk download")
    d = RAW / "gwdb"
    d.mkdir(exist_ok=True)
    z = get(GWDB, d / "GWDBDownload.zip")
    want = ("WaterLevelsMajor.txt", "WellMain.txt")
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            if n.split("/")[-1] in want:
                target = d / n.split("/")[-1]
                if not target.exists():
                    with zf.open(n) as src, open(target, "wb") as dst:
                        dst.write(src.read())
                print(f"  extracted {target.name}")

    (RAW / "download_log.json").write_text(json.dumps(log, indent=1))
    print("done")


if __name__ == "__main__":
    sys.exit(main())
