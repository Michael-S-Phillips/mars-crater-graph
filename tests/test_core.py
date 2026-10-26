"""Sanity tests for geometry, pair classification and edge direction."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.catalog import destination, bearing, frac_circle_inside, R_MARS_KM
from mcgraph.pairs import find_pairs
from mcgraph.graph import edge_table


def test_destination_bearing_roundtrip():
    la, lo = destination(-20.0, 330.0, np.radians(60), 100 / R_MARS_KM)
    b = bearing(-20.0, 330.0, la, lo)
    assert abs(np.degrees(b) - 60) < 1e-6


def test_frac_circle_inside():
    assert frac_circle_inside(1, 10, 0) == 1.0           # small fully inside
    assert frac_circle_inside(1, 10, 20) == 0.0          # disjoint
    # centre of B on A's rim, B much smaller -> about half of B inside
    assert abs(frac_circle_inside(0.01, 10, 10) - 0.5) < 0.01


def _cat(rows):
    df = pd.DataFrame(rows, columns=["lat", "lon", "D_km"])
    df["r_km"] = df.D_km / 2; df["ang_r"] = df.r_km / R_MARS_KM
    la, lo = np.radians(df.lat), np.radians(df.lon)
    df["x"], df["y"], df["z"] = np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)
    return df


def test_pair_kinds_and_containment_direction():
    k = R_MARS_KM * np.pi / 180  # km per degree
    df = _cat([[0, 10, 100],                 # 0: big crater
               [0, 10 + 10 / k, 10],         # 1: inside 0
               [0, 10 + 50 / k, 20],         # 2: straddles 0's rim
               [0, 10 + 90 / k, 10]])        # 3: outside the rim, within one radius of it
    p = find_pairs(df)
    kinds = {(int(r.big), int(r.small)): r.kind for r in p.itertuples()}
    assert kinds[(0, 1)] == "contain"
    assert kinds[(0, 2)] == "overlap"
    assert kinds[(0, 3)] == "ejecta_zone"
    p["vote"] = 0
    e = edge_table(p)
    row = e[(e.method == "containment")].iloc[0]
    assert (row.older, row.younger) == (0, 1)


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("ok", name)
