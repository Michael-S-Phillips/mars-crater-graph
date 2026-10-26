"""Extract THEMIS daytime IR (100 m) and HRSC-MOLA blended DEM (200 m) chips
around crater-pair windows, for visual validation.

Run where the global mosaics live:
  python extract_chips.py windows.json OUTDIR THEMIS_COG DEM_COG
Each window: {pair_id, clat, clon (0-360E), hw_km}. Output: OUTDIR/<pair_id>.npz
with isotropic chips (themis 320x320 uint8, dem 160x160 float32, meters).
"""
import json, os, sys
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling

R_M = 3396190.0
wins = json.load(open(sys.argv[1])); out = sys.argv[2]
themis_p, dem_p = sys.argv[3], sys.argv[4]
os.makedirs(out, exist_ok=True)


def read_chip(ds, clat, clon, hw_km, n, degrees):
    # source x = 0 sits at the CRS central meridian (180 for the USGS THEMIS mosaic)
    cm = float(ds.crs.to_dict().get("lon_0", 0.0))
    lon = ((clon - cm + 180) % 360) - 180
    dlat = np.degrees(hw_km * 1000 / R_M)
    dlon = dlat / max(np.cos(np.radians(clat)), 0.05)
    if degrees:  # dataset in degrees
        b = (lon - dlon, clat - dlat, lon + dlon, clat + dlat)
    else:        # simple cylindrical meters
        k = np.radians(1) * R_M
        b = ((lon - dlon) * k, (clat - dlat) * k, (lon + dlon) * k, (clat + dlat) * k)
    w = from_bounds(*b, transform=ds.transform)
    return ds.read(1, window=w, out_shape=(n, n), resampling=Resampling.bilinear,
                   boundless=True, fill_value=0)


with rasterio.open(themis_p) as th, rasterio.open(dem_p) as dm:
    for w in wins:
        f = os.path.join(out, f"{w['pair_id']}.npz")
        if os.path.exists(f):
            continue
        t = read_chip(th, w["clat"], w["clon"], w["hw_km"], 320, False)
        d = read_chip(dm, w["clat"], w["clon"], w["hw_km"], 160, True).astype(np.float32)
        np.savez_compressed(f, themis=t, dem=d, **{k: np.float64(v) for k, v in w.items()})
print("ok", len(wins))
