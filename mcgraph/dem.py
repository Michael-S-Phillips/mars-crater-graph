"""Global topography sampler.

Default source: MGS MOLA MEGDR 128 pixel/degree topography tiles (megt*hb.img,
PDS Geosciences Node), mosaicked once into a single little-endian int16 memmap
covering 88N-88S, 0-360E (simple cylindrical, ~463 m/px at the equator).

Any other simple-cylindrical global DEM (e.g. the HRSC-MOLA blended 200 m DEM) can
be used through `GeoTiffDEM`, which needs rasterio.
"""
from __future__ import annotations

import os
import numpy as np

MOLA_BASE = ("https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/"
             "mgsl_300x/meg128/")
TILE_LATS = ["88n", "44n", "00n", "44s"]
TILE_LONS = ["000", "090", "180", "270"]
PPD = 128
TILE_LINES, TILE_SAMPLES = 44 * PPD, 90 * PPD


def build_mola_mosaic(tile_dir: str, out_path: str) -> str:
    """Mosaic the 16 MEGDR topography tiles into one memmap file."""
    lines, samples = 4 * TILE_LINES, 4 * TILE_SAMPLES
    mm = np.memmap(out_path, dtype="<i2", mode="w+", shape=(lines, samples))
    for a, la in enumerate(TILE_LATS):
        for b, lo in enumerate(TILE_LONS):
            f = os.path.join(tile_dir, f"megt{la}{lo}hb.img")
            t = np.fromfile(f, dtype=">i2").reshape(TILE_LINES, TILE_SAMPLES)
            mm[a * TILE_LINES:(a + 1) * TILE_LINES,
               b * TILE_SAMPLES:(b + 1) * TILE_SAMPLES] = t
    mm.flush()
    return out_path


class MolaDEM:
    """Bilinear sampler on the mosaicked MEGDR (meters relative to the areoid)."""

    lat_max = 88.0
    km_per_px_equator = 2 * np.pi * 3396.0 / (360 * PPD)

    def __init__(self, mosaic_path: str):
        self.a = np.memmap(mosaic_path, dtype="<i2", mode="r",
                           shape=(4 * TILE_LINES, 4 * TILE_SAMPLES))
        self.nl, self.ns = self.a.shape

    def sample(self, lat, lon):
        lat = np.asarray(lat, float); lon = np.mod(np.asarray(lon, float), 360.0)
        r = (self.lat_max - lat) * PPD - 0.5
        c = lon * PPD - 0.5
        r = np.clip(r, 0, self.nl - 1.001)
        r0 = np.floor(r).astype(np.int64); fr = r - r0
        c0 = np.floor(c).astype(np.int64); fc = c - c0
        c0 = np.mod(c0, self.ns); c1 = np.mod(c0 + 1, self.ns)
        a = self.a
        z00 = a[r0, c0].astype(np.float32); z01 = a[r0, c1].astype(np.float32)
        z10 = a[r0 + 1, c0].astype(np.float32); z11 = a[r0 + 1, c1].astype(np.float32)
        return ((1 - fr) * ((1 - fc) * z00 + fc * z01) + fr * ((1 - fc) * z10 + fc * z11))


class GeoTiffDEM:
    """Sampler for a global simple-cylindrical GeoTIFF in meters (needs rasterio).

    Reads the whole raster into memory on first use, so use an overview or a
    regional subset for very large files.
    """

    def __init__(self, path: str, radius_km: float = 3396.19):
        import rasterio
        with rasterio.open(path) as ds:
            self.a = ds.read(1, out_dtype="float32")
            t = ds.transform
            self.x0, self.dx, self.y0, self.dy = t.c, t.a, t.f, t.e
        self.m_per_deg = 2 * np.pi * radius_km * 1000 / 360.0
        self.km_per_px_equator = abs(self.dx) / 1000.0
        self.nl, self.ns = self.a.shape
        self.lat_max = 90.0

    def sample(self, lat, lon):
        lon = np.asarray(lon, float)
        lon = np.where(lon > 180, lon - 360, lon)
        x = lon * self.m_per_deg; y = np.asarray(lat, float) * self.m_per_deg
        c = (x - self.x0) / self.dx - 0.5
        r = (y - self.y0) / self.dy - 0.5
        r = np.clip(r, 0, self.nl - 1.001)
        r0 = np.floor(r).astype(np.int64); fr = r - r0
        c0 = np.floor(c).astype(np.int64); fc = c - c0
        c0 = np.mod(c0, self.ns); c1 = np.mod(c0 + 1, self.ns)
        a = self.a
        return ((1 - fr) * ((1 - fc) * a[r0, c0] + fc * a[r0, c1])
                + fr * ((1 - fc) * a[r0 + 1, c0] + fc * a[r0 + 1, c1]))
