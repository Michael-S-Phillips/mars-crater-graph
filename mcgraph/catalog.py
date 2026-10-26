"""Crater catalog loading and spherical-geometry helpers.

The default catalog is the Robbins & Hynek (2012) Mars crater database in its
2020 release ("Catalog_Mars_Release_2020_1kmPlus_FullMorphData.csv"), which adds
rim, ejecta and floor degradation states (DEG_RIM, DEG_EJC, DEG_FLR; 4 = freshest,
1 = most degraded) to the original morphology fields.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

R_MARS_KM = 3389.5  # IAU mean radius

KEEP = [
    "CRATER_ID", "LAT_CIRC_IMG", "LON_CIRC_IMG", "DIAM_CIRC_IMG", "DIAM_CIRC_SD_IMG",
    "DIAM_ELLI_ELLIP_IMG", "ARC_IMG", "PTS_RIM_IMG", "LAY_NUMBER", "LAY_MORPH1",
    "INT_MORPH1", "INT_MORPH3", "DEG_RIM", "DEG_EJC", "DEG_FLR",
]


def load_robbins(path: str) -> pd.DataFrame:
    """Load the Robbins catalog and add derived columns.

    Adds: lat, lon (0-360 E), D_km, r_km, ang_r (crater radius in radians),
    x, y, z (unit vectors), layered (bool: mapped layered ejecta).
    """
    df = pd.read_csv(path, encoding="utf-8-sig", low_memory=False)
    df = df[[c for c in KEEP if c in df.columns]].copy()
    df = df.rename(columns={"LAT_CIRC_IMG": "lat", "LON_CIRC_IMG": "lon",
                            "DIAM_CIRC_IMG": "D_km"})
    df["lon"] = np.mod(df["lon"], 360.0)
    df["r_km"] = df["D_km"] / 2.0
    df["ang_r"] = df["r_km"] / R_MARS_KM
    lat = np.radians(df["lat"].to_numpy())
    lon = np.radians(df["lon"].to_numpy())
    df["x"] = np.cos(lat) * np.cos(lon)
    df["y"] = np.cos(lat) * np.sin(lon)
    df["z"] = np.sin(lat)
    lm = df.get("LAY_MORPH1")
    if lm is not None:
        lm = lm.fillna("").astype(str)
        df["layered"] = lm.str.contains("LE", regex=False)
    df = df.reset_index(drop=True)
    df.index.name = "node"
    return df


def unit_vectors(df: pd.DataFrame) -> np.ndarray:
    return df[["x", "y", "z"]].to_numpy()


def central_angle(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Great-circle angle (radians) between rows of unit vectors u and v."""
    cross = np.linalg.norm(np.cross(u, v), axis=-1)
    dot = np.sum(u * v, axis=-1)
    return np.arctan2(cross, dot)


def destination(lat_deg, lon_deg, bearing_rad, ang_dist_rad):
    """Spherical direct problem: point at angular distance along a bearing.

    All inputs broadcast. Returns (lat_deg, lon_deg in 0-360).
    """
    lat1 = np.radians(lat_deg)
    lon1 = np.radians(lon_deg)
    sl, cl = np.sin(lat1), np.cos(lat1)
    sd, cd = np.sin(ang_dist_rad), np.cos(ang_dist_rad)
    lat2 = np.arcsin(np.clip(sl * cd + cl * sd * np.cos(bearing_rad), -1, 1))
    lon2 = lon1 + np.arctan2(np.sin(bearing_rad) * sd * cl, cd - sl * np.sin(lat2))
    return np.degrees(lat2), np.mod(np.degrees(lon2), 360.0)


def bearing(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Initial great-circle bearing (radians, clockwise from north)."""
    lat1, lat2 = np.radians(lat1_deg), np.radians(lat2_deg)
    dlon = np.radians(lon2_deg - lon1_deg)
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    return np.arctan2(x, y)


def frac_circle_inside(b, a, d):
    """Fraction of the circumference of circle B (radius b) lying inside circle A
    (radius a) when their centers are separated by d (planar approximation)."""
    b = np.asarray(b, float); a = np.asarray(a, float); d = np.asarray(d, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = (d ** 2 + b ** 2 - a ** 2) / (2 * d * b)
    c = np.where(d == 0, np.where(b <= a, -1.0, 1.0), c)
    return np.arccos(np.clip(c, -1.0, 1.0)) / np.pi
