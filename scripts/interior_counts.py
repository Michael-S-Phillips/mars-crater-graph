"""Interior crater counts and crater-time exposure for every crater >= 10 km.

  python scripts/interior_counts.py --work WORK --out outputs
Writes outputs/interior_counts.parquet (one row per host crater) and
outputs/interior_members.parquet (each counted crater and its host).
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.ages import interior_counts

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--dmin-host", type=float, default=10.0)
ap.add_argument("--dmax-host", type=float, default=600.0)
a = ap.parse_args()
cat = pd.read_parquet(os.path.join(a.work, "catalog.parquet"))
pairs = pd.read_parquet(os.path.join(a.out, "pairs_evidence.parquet"),
                        columns=["big", "small", "near_concentric", "possible_ring"])
e = pd.read_csv(os.path.join(a.out, "edges.csv"))
e = e[e.kept]
hosts = np.where((cat.D_km >= a.dmin_host) & (cat.D_km <= a.dmax_host))[0]
t0 = time.time()
mem = []
r = interior_counts(cat, pairs, e, hosts, dmins=(1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0), members=mem)
r.insert(0, "CRATER_ID", cat.CRATER_ID.to_numpy()[r.index])
r.to_parquet(os.path.join(a.out, "interior_counts.parquet"))
pd.DataFrame(mem, columns=["host", "crater", "D_km"]).to_parquet(
    os.path.join(a.out, "interior_members.parquet"))
print(len(r), "hosts in", round(time.time() - t0), "s")
