"""Crater-count ages for individual craters, measured on their own interiors.

Crater time
-----------
Ages are first expressed as "crater time": tau = N(>=1 km) per km^2 accumulated
since a surface formed. With a spatially uniform impact flux, craters arrive at a
constant rate in tau, whatever the impact-rate history, so counting statistics and
superposition order can be combined without choosing a chronology. A chronology
function converts tau to years only at the end.

Production and chronology functions (coefficients as distributed with craterstats,
Michael et al.):
  production  Ivanov (2001), polynomial in log10 D, valid 0.015-362 km
  chronology  Hartmann & Neukum (2001):  N1(T) = 3.22e-14 (exp(6.93 T) - 1) + 4.875e-4 T
              Ivanov (2001):             N1(T) = 2.68e-14 (exp(6.93 T) - 1) + 4.13e-4 T
  epochs      Michael (2013) densities for the Tanaka (1986) epoch boundaries

Interior counts
---------------
For a host crater H, the count area is H's interior (the disk inside the rim).
Craters whose centers fall inside H postdate it. A crater of diameter D can only
be recorded where no larger younger crater has since excavated the surface, so
for each diameter the area excludes the interiors of larger superposed craters,
and a crater is counted only if it is not itself inside a larger superposed crater
(the non-sparseness correction of Kneissl et al., 2016). With the production
function's shape fixed, the count k of craters >= Dmin is sufficient for tau, and

    k ~ Poisson(tau * Lambda),   Lambda = integral_Dmin^D_H  f(D) A_eff(D) dD,

where f is the production function's differential shape normalized to N(>=1) = 1.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .catalog import R_MARS_KM

PF_IVANOV = [-3.383677e00, -3.197453e00, 1.256814e00, 7.915374e-01, -4.860814e-01,
             -3.630098e-01, 1.015683e-01, 6.755923e-02, -1.180639e-02, -4.753462e-03,
             6.232845e-04, 5.805492e-05]
CF = {"Hartmann & Neukum (2001)": (3.22e-14, 6.93, 4.875e-4),
      "Ivanov (2001)": (2.68e-14, 6.93, 4.13e-4)}
DEFAULT_CF = "Hartmann & Neukum (2001)"

# Michael (2013) epoch boundaries: cumulative density N(>=D_ref) per 10^6 km^2
EPOCHS = [("LA/MA", 1, 160), ("MA/EA", 1, 600), ("EA/LH", 1, 2100), ("LH/EH", 5, 125),
          ("EH/LN", 5, 200), ("LN/MN", 16, 100), ("MN/EN", 16, 200)]


def pf_rel(D):
    """N(>=D) / N(>=1 km) for the Ivanov (2001) production function."""
    x = np.log10(np.asarray(D, dtype=float))
    s = sum(a * x ** n for n, a in enumerate(PF_IVANOV) if n > 0)
    return 10.0 ** s


def n1_of_t(T, cf=DEFAULT_CF):
    a, b, c = CF[cf]
    return a * (np.exp(b * np.asarray(T, dtype=float)) - 1.0) + c * np.asarray(T, dtype=float)


_TGRID = np.linspace(0.0, 4.6, 46001)


def t_of_n1(n1, cf=DEFAULT_CF):
    """Model age (Ga) for crater time tau = N(>=1)/km^2."""
    return np.interp(np.asarray(n1, dtype=float), n1_of_t(_TGRID, cf), _TGRID)


def epoch_boundaries_tau():
    """Epoch boundaries in crater time (N1 per km^2)."""
    return {name: dens / 1e6 / pf_rel(dref) for name, dref, dens in EPOCHS}


def _local_xy(lat0, lon0, lat, lon):
    """Azimuthal-equidistant coordinates (km) about (lat0, lon0)."""
    p0, l0 = np.radians(lat0), np.radians(lon0)
    p, l = np.radians(lat), np.radians(lon)
    cosc = np.sin(p0) * np.sin(p) + np.cos(p0) * np.cos(p) * np.cos(l - l0)
    c = np.arccos(np.clip(cosc, -1, 1))
    az = np.arctan2(np.sin(l - l0) * np.cos(p),
                    np.cos(p0) * np.sin(p) - np.sin(p0) * np.cos(p) * np.cos(l - l0))
    d = c * R_MARS_KM
    return d * np.sin(az), d * np.cos(az)


def _disk_points(n=1500, seed=0):
    """Fixed quasi-uniform points in the unit disk (Fermat spiral)."""
    i = np.arange(n) + 0.5
    r = np.sqrt(i / n)
    a = i * np.pi * (3 - np.sqrt(5))
    return r * np.cos(a), r * np.sin(a)


def interior_counts(cat: pd.DataFrame, pairs: pd.DataFrame, edges: pd.DataFrame,
                    hosts, dmins=(1.0, 1.5, 2.0, 3.0), n_pts=1500, dgrid_n=60, members=None):
    """Interior crater counts and exposure Lambda for each host crater.

    cat    catalog with lat, lon, D_km, x, y, z (index = node)
    pairs  pair table (big, small, near_concentric, possible_ring)
    edges  kept graph edges (older, younger)
    hosts  node indices of host craters
    Returns a DataFrame indexed by host node with k_<Dmin>, L_<Dmin> (km^2 per unit
    tau), the host area, and the number of excluded superposed craters.
    members  optional list; if given, (host, crater, D_km) is appended for every
             counted crater >= min(dmins)
    """
    xyz = cat[["x", "y", "z"]].to_numpy()
    lat, lon, D = cat.lat.to_numpy(), cat.lon.to_numpy(), cat.D_km.to_numpy()
    tree = cKDTree(xyz)
    flagged = set(zip(pairs.big[pairs.near_concentric | pairs.possible_ring],
                      pairs.small[pairs.near_concentric | pairs.possible_ring]))
    older_of = {}   # host -> set of craters the graph says are OLDER than it
    younger_of = {}  # host -> set of craters the graph says are younger
    for o, y in zip(edges.older.to_numpy(), edges.younger.to_numpy()):
        younger_of.setdefault(o, set()).add(y)
        older_of.setdefault(y, set()).add(o)
    px, py = _disk_points(n_pts)
    rows = []
    for h in hosts:
        rh = D[h] / 2
        theta = rh / R_MARS_KM
        area = 2 * np.pi * R_MARS_KM ** 2 * (1 - np.cos(theta))
        # candidates out to 1.5 radii, so overlapping younger craters are found
        cand = np.array(tree.query_ball_point(xyz[h], 2 * np.sin(min(np.pi, 1.5 * theta + 0.02) / 2)))
        cand = cand[(cand != h) & (D[cand] < D[h])]
        x, y = _local_xy(lat[h], lon[h], lat[cand], lon[cand])
        dist = np.hypot(x, y)
        old = older_of.get(h, set())
        yng = younger_of.get(h, set())
        bad = np.array([(h, j) in flagged or j in old for j in cand], dtype=bool)
        inside = (dist < rh) & ~bad
        outside_young = (dist >= rh) & (dist < rh + D[cand] / 2) & np.array(
            [j in yng for j in cand], dtype=bool)
        excl = inside | outside_young
        ex_x, ex_y, ex_r, ex_D = x[excl], y[excl], D[cand][excl] / 2, D[cand][excl]
        # a counted crater must not lie inside a larger superposed crater
        ix, iy, iD = x[inside], y[inside], D[cand][inside]
        if len(ix):
            dd = np.hypot(ix[:, None] - ex_x[None, :], iy[:, None] - ex_y[None, :])
            covered = ((dd < ex_r[None, :]) & (ex_D[None, :] > iD[:, None])).any(axis=1)
        else:
            covered = np.zeros(0, bool)
        # A_eff(D): fraction of the interior not inside a superposed crater larger than D
        big = ex_D >= min(dmins)
        qx, qy = px * rh, py * rh
        if big.any():
            inn = (np.hypot(qx[:, None] - ex_x[big][None, :], qy[:, None] - ex_y[big][None, :])
                   < ex_r[big][None, :])
            dcover = np.where(inn, ex_D[big][None, :], 0.0).max(axis=1)
        else:
            dcover = np.zeros(len(qx))
        dc_sorted = np.sort(dcover)
        row = {"node": h, "D_km": D[h], "area_km2": area, "n_inside": int(inside.sum()),
               "n_excl_ge1": int(big.sum())}
        for dm in dmins:
            if dm >= D[h]:
                row[f"k_{dm:g}"], row[f"L_{dm:g}"] = 0, 0.0
                continue
            edges_d = np.geomspace(dm, D[h], dgrid_n + 1)
            mid = np.sqrt(edges_d[:-1] * edges_d[1:])
            frac = np.searchsorted(dc_sorted, mid, side="right") / len(dc_sorted)
            w = pf_rel(edges_d[:-1]) - pf_rel(edges_d[1:])
            row[f"L_{dm:g}"] = float(area * np.sum(w * frac))
            row[f"k_{dm:g}"] = int(((iD >= dm) & ~covered).sum())
        rows.append(row)
        if members is not None:
            cnt = (iD >= min(dmins)) & ~covered
            members.extend(zip([h] * int(cnt.sum()), cand[inside][cnt], iD[cnt]))
    return pd.DataFrame(rows).set_index("node")


def gamma_post(k, L):
    """Posterior of tau under a uniform prior: Gamma(k + 1, rate L). Returns shape, rate."""
    return np.asarray(k, float) + 1.0, np.asarray(L, float)
