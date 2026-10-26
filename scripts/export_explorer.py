"""Export compact JSON + a shaded-relief image for the interactive explorer page.

  python scripts/export_explorer.py --work WORK --out outputs --dest explorer
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.dem import MolaDEM, PPD

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--dest", required=True)
a = ap.parse_args()
os.makedirs(os.path.join(a.dest, "data"), exist_ok=True)

n = pd.read_csv(os.path.join(a.out, "nodes.csv.gz"), index_col="node", low_memory=False)
e = pd.read_csv(os.path.join(a.out, "edges.csv"))
g = n[n.in_graph].copy()
idx = {node: i for i, node in enumerate(g.index)}
units = sorted(g.unit.dropna().unique().tolist())
uidx = {u: i for i, u in enumerate(units)}
nodes = {
    "id": g.CRATER_ID.tolist(),
    "lat": (g.lat * 1e4).round().astype(int).tolist(),
    "lon": (g.lon * 1e4).round().astype(int).tolist(),
    "D": (g.D_km * 100).round().astype(int).tolist(),
    "depth": g.depth.astype(int).tolist(),
    "height": g.height.astype(int).tolist(),
    "nOld": g.n_older.astype(int).tolist(),
    "nYoung": g.n_younger.astype(int).tolist(),
    "deg": g.DEG_RIM.fillna(0).astype(int).tolist(),
    "unit": [uidx.get(u, -1) if isinstance(u, str) else -1 for u in g.unit],
}
desc = g.drop_duplicates("unit").set_index("unit").unit_desc.to_dict()
meth = {"containment": 0, "arc": 1, "topo+arc": 2, "topo": 3}
h = e[e.kept & e.in_hasse]
edges = {"o": [idx[x] for x in h.older], "y": [idx[x] for x in h.younger],
         "m": [meth[m] for m in h.method]}
summary = json.load(open(os.path.join(a.out, "summary.json")))
val = json.load(open(os.path.join(a.out, "validation.json")))
meta = {"units": units, "unitDesc": [desc.get(u, "") for u in units],
        "methods": ["containment", "rim-arc", "rim-arc + topography", "topography"],
        "summary": summary, "validation": val}
# craters with no relation in the default graph, and why (the most specific reason)
#   0 no other cataloged crater touches it        1 only on or near another crater's ejecta
#   2 rims barely touch another crater           3 nested pair set aside (duplicate or basin ring)
#   4 rims cross, but neither test decided       5 only a topography-only (low-confidence) relation
p = pd.read_parquet(os.path.join(a.out, "pairs_evidence.parquet"))
code = pd.Series(0, index=n.index, dtype=int)
for mask, val in (
        (p.kind == "ejecta_zone", 1), (p.kind == "touch", 2),
        ((p.kind == "contain") & (p.near_concentric | p.possible_ring), 3),
        ((p.kind == "overlap") & (p.vote.fillna(0) == 0), 4),
        ((p.kind == "overlap") & (p.confidence == "low"), 5)):
    hit = np.union1d(p.big[mask], p.small[mask])
    code.loc[hit] = np.maximum(code.loc[hit], val)
u = n[~n.in_graph]
unlinked = {
    "id": u.CRATER_ID.tolist(),
    "lat": (u.lat * 1e3).round().astype(int).tolist(),
    "lon": (u.lon * 1e3).round().astype(int).tolist(),
    "D": (u.D_km * 100).round().astype(int).tolist(),
    "st": code.loc[u.index].astype(int).tolist(),
}
json.dump(unlinked, open(os.path.join(a.dest, "data", "unlinked.json"), "w"), separators=(",", ":"))
json.dump(nodes, open(os.path.join(a.dest, "data", "nodes.json"), "w"), separators=(",", ":"))
json.dump(edges, open(os.path.join(a.dest, "data", "edges.json"), "w"), separators=(",", ":"))
json.dump(meta, open(os.path.join(a.dest, "data", "meta.json"), "w"), separators=(",", ":"))

# shaded relief, 8 px/degree, 88N-88S, 0-360E
dem = MolaDEM(os.path.join(a.work, "mola128_global.i2"))
step = 16
z = np.asarray(dem.a[::step, ::step], dtype=float)
px = 3396000 * np.pi / 180 / (PPD / step)
lat = np.linspace(88, -88, z.shape[0])
gy, gx = np.gradient(z, px, 1.0)
gx = gx / (px * np.cos(np.radians(lat))[:, None])
slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
A, E = np.radians(315), np.radians(40)
hs = np.sin(E) * np.cos(slope) + np.cos(E) * np.sin(slope) * np.cos(A - aspect)
zn = (z - np.percentile(z, 1)) / (np.percentile(z, 99) - np.percentile(z, 1))
img = np.clip(0.75 * hs + 0.25 * np.clip(zn, 0, 1), 0, 1)
from PIL import Image
Image.fromarray((img * 255).astype(np.uint8)).save(os.path.join(a.dest, "data", "relief.jpg"), quality=82)
for f in os.listdir(os.path.join(a.dest, "data")):
    print(f, os.path.getsize(os.path.join(a.dest, "data", f)))
