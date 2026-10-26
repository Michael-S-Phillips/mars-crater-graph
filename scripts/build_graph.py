"""Step 3: build the directed stratigraphic graph and export it.

  python scripts/build_graph.py --work WORKDIR --geology SIM3292_Global_Geology.shp --out outputs
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.graph import edge_table, build_graph, break_cycles, node_metrics, hasse_edges
from mcgraph.units import host_units

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True)
ap.add_argument("--geology", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-conf", default="medium", choices=["high", "medium", "low"])
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()

df = pd.read_parquet(os.path.join(a.work, "catalog.parquet"))
pairs = pd.read_parquet(os.path.join(a.work, "pairs_evidence.parquet"))
units = host_units(df, a.geology)
print(f"units joined ({time.time()-t0:.0f}s)")

E = edge_table(pairs, a.min_conf)
G = build_graph(E)
n_edges_raw = G.number_of_edges()
removed = break_cycles(G)
print(f"edges {n_edges_raw}, removed for cycles {len(removed)} ({time.time()-t0:.0f}s)")
assert nx.is_directed_acyclic_graph(G)
M = node_metrics(G)
R = hasse_edges(G)
print(f"metrics + transitive reduction ({time.time()-t0:.0f}s); hasse edges {R.number_of_edges()}")

cid = df["CRATER_ID"]
# ---- node table
nodes = df[["CRATER_ID", "lat", "lon", "D_km", "ARC_IMG", "DEG_RIM", "DEG_EJC", "DEG_FLR",
            "LAY_MORPH1"]].join(units).join(M, how="left")
nodes["in_graph"] = nodes["depth"].notna()
for c in ["depth", "height", "n_older", "n_younger", "component", "component_size"]:
    nodes[c] = nodes[c].astype("Int64")
nodes.to_csv(os.path.join(a.out, "nodes.csv.gz"), index_label="node")

# ---- edge tables (CRATER_IDs alongside node indices)
alive = {(u, v) for u, v in G.edges}
E["kept"] = [(u, v) in alive for u, v in zip(E.older, E.younger)]
E["in_hasse"] = [R.has_edge(u, v) for u, v in zip(E.older, E.younger)]
E.insert(1, "older_id", cid.to_numpy()[E.older]); E.insert(2, "younger_id", cid.to_numpy()[E.younger])
E.to_csv(os.path.join(a.out, "edges.csv"), index=False)
pd.DataFrame(removed, columns=["older", "younger", "cycle_len", "pair_id", "kind", "method", "confidence"]).to_csv(
    os.path.join(a.out, "edges_removed_cycles.csv"), index=False)

# ---- unresolved overlaps (candidates for imagery-based classification)
un = pairs[(pairs.kind == "overlap") & (pairs.vote.fillna(0) == 0)].copy()
un.insert(0, "big_id", cid.to_numpy()[un.big]); un.insert(1, "small_id", cid.to_numpy()[un.small])
un[["big_id", "small_id", "big", "small", "d_km", "r_big", "r_small", "confidence",
    "P_big", "P_small", "y_big", "y_small"]].to_csv(
    os.path.join(a.out, "unresolved_overlaps.csv"), index_label="pair_id")
pairs.to_parquet(os.path.join(a.out, "pairs_evidence.parquet"))

# ---- GraphML (full DAG, crater attributes on nodes)
H = nx.DiGraph()
attrs = nodes.loc[list(G.nodes)]
for n, r in attrs.iterrows():
    H.add_node(r.CRATER_ID, lat=float(r.lat), lon=float(r.lon), D_km=float(r.D_km),
               depth=int(r.depth), height=int(r.height), unit=str(r.unit))
for u, v, d in G.edges(data=True):
    H.add_edge(cid.iat[u], cid.iat[v], kind=d["kind"], method=d["method"],
               confidence=d["confidence"], hasse=R.has_edge(u, v))
nx.write_graphml(H, os.path.join(a.out, "mars_crater_strat_graph.graphml.gz"))

# ---- summary
comp_sizes = nodes.drop_duplicates("component").component_size.dropna().astype(int)
summary = {
    "craters_in_catalog": int(len(df)),
    "pairs": pairs.kind.value_counts().to_dict(),
    "overlap_confidence": pairs[pairs.kind == "overlap"].confidence.value_counts().to_dict(),
    "edges_by_method": E[E.kept].method.value_counts().to_dict(),
    "edges_total": int(E.kept.sum()),
    "edges_removed_for_cycles": len(removed),
    "hasse_edges": R.number_of_edges(),
    "craters_in_graph": int(nodes.in_graph.sum()),
    "components": int(len(comp_sizes)),
    "largest_components": comp_sizes.sort_values(ascending=False).head(10).tolist(),
    "max_depth": int(nodes.depth.max()),
    "depth_counts": nodes.depth.value_counts().sort_index().astype(int).to_dict(),
}
summary["depth_counts"] = {int(k): int(v) for k, v in summary["depth_counts"].items()}
json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=2)
print(json.dumps(summary, indent=2))
print(f"done ({time.time()-t0:.0f}s)")
