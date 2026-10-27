"""Cut a lat/lon box out of the global THEMIS daytime IR mosaic, resampled to a
regular grid in degrees (rows from north, columns from west), for figure basemaps.

  python scripts/extract_themis_cutout.py SRC.tif OUT.npz LAT0 LAT1 LON0 LON1 --ppd 512

LON0 > LON1 means the box crosses 0 E. Run where the mosaic lives. The USGS mosaic is
equirectangular in meters with its central meridian at 180 E; the script reads the
central meridian from the CRS, as make_tiles.py does.
"""
import argparse, math
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("lat0", type=float); ap.add_argument("lat1", type=float)
ap.add_argument("lon0", type=float); ap.add_argument("lon1", type=float)
ap.add_argument("--ppd", type=float, default=512)
a = ap.parse_args()
lat0, lat1 = min(a.lat0, a.lat1), max(a.lat0, a.lat1)
lon0, lon1 = a.lon0 % 360, a.lon1 % 360
if lon1 <= lon0:
    lon1 += 360
with rasterio.open(a.src) as ds:
    k = math.radians(1) * ds.crs.to_dict().get("R", ds.crs.to_dict().get("a", 3396190.0))
    cm = float(ds.crs.to_dict().get("lon_0", 0.0))
    h = int(round((lat1 - lat0) * a.ppd))
    parts, lo = [], lon0
    while lo < lon1 - 1e-9:
        x_w = ((lo - cm + 180) % 360) - 180          # source x (degrees) of this piece's west edge
        span = min(lon1 - lo, 180 - x_w)             # stop at the source grid's east edge
        w = from_bounds(x_w * k, lat0 * k, (x_w + span) * k, lat1 * k, transform=ds.transform)
        parts.append(ds.read(1, window=w, out_shape=(h, int(round(span * a.ppd))),
                             resampling=Resampling.average))
        lo += span
    img = np.concatenate(parts, axis=1)
np.savez_compressed(a.out, img=img, lat0=lat0, lat1=lat1, lon0=lon0, lon1=lon1, ppd=a.ppd)
print(a.out, img.shape, img.dtype, int(img.min()), int(img.max()))
