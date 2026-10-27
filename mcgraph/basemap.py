"""Figure basemaps: THEMIS daytime IR blended with MOLA shaded relief.

THEMIS cutouts come from scripts/extract_themis_cutout.py (a regular lat/lon grid,
rows from north). The blend is a weighted mean of the two images after a
percentile stretch, so the 100 m THEMIS texture carries the detail and the 463 m
MOLA hillshade carries the relief; both stay grayscale under colored overlays.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates


def load_cutout(path):
    z = np.load(path)
    return {k: z[k] for k in z.files}


def sample_cutout(c, lat, lon):
    """Bilinear sample of a cutout at lat/lon arrays (degrees, any longitude convention)."""
    ppd = float(c["ppd"])
    row = (float(c["lat1"]) - np.asarray(lat, float)) * ppd - 0.5
    col = ((np.asarray(lon, float) - float(c["lon0"])) % 360) * ppd - 0.5
    img = c["img"].astype(float)
    img[img == 0] = np.nan  # no-data
    out = map_coordinates(np.nan_to_num(img, nan=np.nanmedian(img)), [row.ravel(), col.ravel()],
                          order=1, mode="nearest")
    return out.reshape(np.shape(lat))


def stretch(a, lo=1, hi=99):
    p0, p1 = np.nanpercentile(a, [lo, hi])
    return np.clip((a - p0) / (p1 - p0), 0, 1)


def hillshade(z, dx, dy, az=315, el=40):
    """Hillshade of an elevation grid with rows from north (dx, dy in meters)."""
    gy, gx = np.gradient(z, dy, dx)
    gy = -gy
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
    A, E = np.radians(az), np.radians(el)
    return np.sin(E) * np.cos(slope) + np.cos(E) * np.sin(slope) * np.cos(A - aspect)


def blend(themis, hs, w_themis=0.6):
    """Grayscale blend of a THEMIS image and a hillshade on the same grid, in 0-1."""
    out = w_themis * stretch(themis) + (1 - w_themis) * stretch(hs, 2, 99.5)
    return stretch(out, 0.5, 99.5)
