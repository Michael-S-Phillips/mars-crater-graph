"""Find spatially related crater pairs and classify their geometry.

For every crater i we search neighbors j with r_j <= r_i inside 3 r_i (great-circle
distance), then classify each pair by how the two rim circles relate:

  contain      j lies entirely inside i's rim (d + r_j <= r_i). i must be older:
               i's excavation would have erased any earlier crater inside its rim.
  overlap      the rims cross, with j penetrating i by >= 10 % of r_j. One crater
               truncates the other, but geometry alone cannot say which.
  touch        rims cross by < 10 % of r_j (too marginal to orient).
  ejecta_zone  j sits outside i's rim but within one i-radius of it, i.e. on i's
               continuous ejecta. Orientation needs imagery (not attempted yet).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .catalog import R_MARS_KM, central_angle, frac_circle_inside, unit_vectors


def find_pairs(df: pd.DataFrame, search_factor: float = 3.0,
               min_penetration: float = 0.10) -> pd.DataFrame:
    u = unit_vectors(df)
    ang_r = df["ang_r"].to_numpy()
    r = df["r_km"].to_numpy()
    tree = cKDTree(u)
    # chord length for angular search radius
    search = np.minimum(search_factor * ang_r, np.pi)
    chord = 2 * np.sin(search / 2)
    # process in size order so huge basins don't dominate memory at once
    I, J = [], []
    for start in range(0, len(df), 50_000):
        idx = np.arange(start, min(start + 50_000, len(df)))
        hits = tree.query_ball_point(u[idx], chord[idx])
        for i, hj in zip(idx, hits):
            if len(hj) <= 1:
                continue
            hj = np.asarray(hj)
            # keep j strictly "smaller" than i (ties broken by index)
            keep = (r[hj] < r[i]) | ((r[hj] == r[i]) & (hj > i))
            hj = hj[keep]
            if hj.size:
                I.append(np.full(hj.size, i)); J.append(hj)
    I = np.concatenate(I); J = np.concatenate(J)
    d = central_angle(u[I], u[J]) * R_MARS_KM
    ri, rj = r[I], r[J]
    kind = np.full(I.size, "", dtype=object)
    contain = d + rj <= ri
    cross = (~contain) & (d < ri + rj)
    pen = (ri + rj - d) / rj
    kind[contain] = "contain"
    kind[cross & (pen >= min_penetration)] = "overlap"
    kind[cross & (pen < min_penetration)] = "touch"
    ej = (~contain) & (~cross) & (d < 2 * ri + rj)
    kind[ej] = "ejecta_zone"
    out = pd.DataFrame({
        "big": I, "small": J, "d_km": d, "r_big": ri, "r_small": rj,
        "size_ratio": rj / ri, "kind": kind,
    })
    out = out[out.kind != ""].reset_index(drop=True)
    # fraction of each rim lying inside the other circle
    out["f_big_in_small"] = frac_circle_inside(out.r_big, out.r_small, out.d_km)
    out["f_small_in_big"] = frac_circle_inside(out.r_small, out.r_big, out.d_km)
    out["small_center_inside_big"] = out.d_km < out.r_big
    # near-concentric, near-equal circles are often catalog artifacts (rings,
    # concentric fill, duplicates); flag rather than trust
    out["near_concentric"] = (out.d_km < 0.15 * out.r_big) & (out.size_ratio > 0.7)
    # large, nearly concentric nested circles with ring-like diameter ratios are
    # more likely the peak ring / inner ring of one basin than two craters
    out["possible_ring"] = ((out.kind == "contain") & (out.d_km < 0.2 * out.r_big)
                            & (out.size_ratio >= 0.3) & (out.size_ratio <= 0.9)
                            & (out.r_small >= 10))
    return out
