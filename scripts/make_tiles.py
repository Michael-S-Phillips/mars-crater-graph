"""Cut a global simple-cylindrical raster into a basemap tile pyramid for the explorer.

Tile scheme (equirectangular, east longitude 0-360):
  level z has 8 * 2**z pixels per degree; tiles are 360 x 360 px and span 45 / 2**z degrees
  rows r = 0 .. 4 * 2**z - 1 from 90N southward; columns c = 0 .. 8 * 2**z - 1 from 0E eastward
  file: <out>/<z>/<r>/<c>.jpg   (tiles with no data are not written)

  python make_tiles.py SRC OUT --units meters|degrees --zmin 0 --zmax 4 \
      [--stretch LO HI] [--rows R0 R1] [--quality 80]

Run it where the source mosaic lives. It is resumable: a row whose files exist is skipped.
"""
import argparse, math, os, sys, time  # noqa: E401
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from PIL import Image

TILE = 360
R_M = 3396190.0

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--units", choices=["meters", "degrees"], required=True)
ap.add_argument("--zmin", type=int, default=0); ap.add_argument("--zmax", type=int, default=4)
ap.add_argument("--stretch", type=float, nargs=2, default=None)
ap.add_argument("--rows", type=int, nargs=2, default=None, help="row range [r0, r1) for the top level only")
ap.add_argument("--quality", type=int, default=80)
ap.add_argument("--budget", type=float, default=160, help="stop after this many seconds")
a = ap.parse_args()
t0 = time.time()
k = math.radians(1) * R_M if a.units == "meters" else 1.0


def central_meridian(ds):
    """Longitude at x = 0 of the source grid. USGS THEMIS mosaics use 180,
    the HRSC-MOLA blends and Murray Lab CTX use 0; reading it from the CRS avoids
    a 180-degree shift."""
    try:
        return float(ds.crs.to_dict().get("lon_0", 0.0))
    except Exception:
        return 0.0


def read_strip(ds, lat_top, lat_bot, ppd):
    """Strip covering 0-360E at ppd pixels/degree between two latitudes."""
    h = int(round((lat_top - lat_bot) * ppd))
    cm = central_meridian(ds)
    halves = []
    for lon0 in (0, 180):  # 0..180E, then 180..360E
        x0 = ((lon0 - cm + 180) % 360) - 180  # source x (degrees) of the half's west edge
        w = from_bounds(x0 * k, lat_bot * k, (x0 + 180) * k, lat_top * k, transform=ds.transform)
        halves.append(ds.read(1, window=w, out_shape=(h, int(180 * ppd)),
                              resampling=Resampling.average))  # boundless reads ignore resampling
    return np.concatenate(halves, axis=1)


with rasterio.open(a.src) as ds:
    for z in range(a.zmin, a.zmax + 1):
        ppd, span = 8 * 2 ** z, 45 / 2 ** z
        nrow, ncol = 4 * 2 ** z, 8 * 2 ** z
        rows = range(*a.rows) if (a.rows and z == a.zmax) else range(nrow)
        for r in rows:
            d = os.path.join(a.out, str(z), str(r))
            if os.path.exists(os.path.join(d, ".done")):
                continue
            if time.time() - t0 > a.budget:
                print(f"budget reached at z={z} r={r}"); sys.exit(0)
            lat_top = 90 - r * span
            strip = read_strip(ds, lat_top, lat_top - span, ppd)
            valid = strip > 0
            if a.stretch:
                lo, hi = a.stretch
                strip = np.clip((strip.astype(np.float32) - lo) / (hi - lo) * 255, 1, 255).astype(np.uint8)
            os.makedirs(d, exist_ok=True)
            for c in range(ncol):
                if not valid[:, c * TILE:(c + 1) * TILE].any():
                    continue
                Image.fromarray(strip[:, c * TILE:(c + 1) * TILE]).save(
                    os.path.join(d, f"{c}.jpg"), "JPEG", quality=a.quality, optimize=True)
            open(os.path.join(d, ".done"), "w").close()
        print(f"z={z} done ({time.time() - t0:.0f}s)")
print("finished")
