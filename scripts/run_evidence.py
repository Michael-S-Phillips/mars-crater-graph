"""Step 1-2: find crater pairs and orient overlapping ones.

Usage:
  python scripts/run_evidence.py --catalog CATALOG.csv --mola-tiles DIR --work DIR
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.catalog import load_robbins
from mcgraph.pairs import find_pairs
from mcgraph.dem import build_mola_mosaic, MolaDEM
from mcgraph.evidence import rim_profiles, topo_test, arc_test, deg_compare, combine

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--mola-tiles", required=True, help="dir with megt*hb.img tiles")
ap.add_argument("--work", required=True)
a = ap.parse_args()
os.makedirs(a.work, exist_ok=True)

t0 = time.time()
df = load_robbins(a.catalog)
df.to_parquet(os.path.join(a.work, "catalog.parquet"))
pairs = find_pairs(df)
print(f"pairs: {pairs.kind.value_counts().to_dict()}  ({time.time()-t0:.0f}s)")

mosaic = os.path.join(a.work, "mola128_global.i2")
if not os.path.exists(mosaic):
    build_mola_mosaic(a.mola_tiles, mosaic)
dem = MolaDEM(mosaic)

ov = pairs[pairs.kind == "overlap"].copy()
nodes = np.union1d(ov.big, ov.small)
prof = rim_profiles(df, nodes, dem)
print(f"rim profiles: {len(prof)}  ({time.time()-t0:.0f}s)")
t = topo_test(ov, df, prof)
ar = arc_test(ov, df)
ov = pd.concat([ov, t, ar], axis=1)
ov["deg_vote"] = deg_compare(ov, df)
ov = pd.concat([ov, combine(ov)], axis=1)
pairs = pairs.join(ov[[c for c in ov.columns if c not in pairs.columns]])
pairs.to_parquet(os.path.join(a.work, "pairs_evidence.parquet"))
print(ov.confidence.value_counts().to_dict())
print(f"done ({time.time()-t0:.0f}s)")
