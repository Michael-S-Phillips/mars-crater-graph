"""Assemble the crater superposition graph and compute stratigraphic metrics.

Edges point from OLDER to YOUNGER (u -> v means v formed after u).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx

CONF_RANK = {"high": 3, "medium": 2, "low": 1}


def edge_table(pairs: pd.DataFrame, min_conf: str = "medium") -> pd.DataFrame:
    """Turn classified pairs into a directed edge list (older -> younger)."""
    keep_rank = CONF_RANK[min_conf]
    c = pairs[pairs.kind == "contain"].copy()
    # near-identical concentric circles are usually duplicate catalog entries
    dup = c.near_concentric & (c.size_ratio > 0.9)
    c = c[~dup]
    c["older"], c["younger"] = c.big, c.small
    c["method"] = "containment"
    c["confidence"] = np.where(c.possible_ring, "low",
                               np.where(c.near_concentric, "medium", "high"))
    o = pairs[(pairs.kind == "overlap") & (pairs.vote.fillna(0) != 0)].copy()
    o["older"] = np.where(o.vote > 0, o.big, o.small)
    o["younger"] = np.where(o.vote > 0, o.small, o.big)
    e = pd.concat([c, o], ignore_index=False)
    e["pair_id"] = e.index
    e = e[e.confidence.map(CONF_RANK).fillna(0) >= keep_rank]
    cols = ["pair_id", "older", "younger", "kind", "method", "confidence", "d_km",
            "r_big", "r_small", "size_ratio"]
    return e[cols].reset_index(drop=True)


def build_graph(edges: pd.DataFrame) -> nx.DiGraph:
    G = nx.DiGraph()
    for r in edges.itertuples(index=False):
        G.add_edge(int(r.older), int(r.younger), pair_id=int(r.pair_id), kind=r.kind,
                   method=r.method, confidence=r.confidence)
    return G


def break_cycles(G: nx.DiGraph) -> list[dict]:
    """Remove the weakest edge in each cycle until the graph is acyclic.

    Containment edges are geometric certainties, so overlap edges are removed
    first; ties go to the edge whose two craters are most similar in size.
    """
    removed = []
    while True:
        sccs = [s for s in nx.strongly_connected_components(G) if len(s) > 1]
        if not sccs:
            return removed
        for s in sccs:
            H = G.subgraph(s)
            try:
                cyc = nx.find_cycle(H)
            except nx.NetworkXNoCycle:
                continue
            def weight(e):
                d = G.edges[e]
                return (d["kind"] == "contain", CONF_RANK.get(d["confidence"], 0))
            u, v = min(((u, v) for u, v, *_ in cyc), key=weight)
            removed.append(dict(older=u, younger=v, cycle_len=len(cyc), **G.edges[u, v]))
            G.remove_edge(u, v)


def node_metrics(G: nx.DiGraph) -> pd.DataFrame:
    """Per-crater stratigraphic metrics on a DAG.

    depth   longest chain of craters this crater is superposed on (0 = no known
            older crater beneath it)
    height  longest chain of craters superposed on it (0 = nothing known younger)
    n_older / n_younger   number of craters known (transitively) to be older /
            younger
    component, component_size   weakly connected component id and size
    """
    order = list(nx.topological_sort(G))
    depth = {}
    for n in order:
        depth[n] = max((depth[p] + 1 for p in G.predecessors(n)), default=0)
    height = {}
    for n in reversed(order):
        height[n] = max((height[s] + 1 for s in G.successors(n)), default=0)
    comp, csize = {}, {}
    for i, c in enumerate(sorted(nx.weakly_connected_components(G), key=len, reverse=True)):
        for n in c:
            comp[n] = i; csize[n] = len(c)
    # transitive counts via bitsets per component (components are modest in size)
    n_older, n_younger = {}, {}
    for c in nx.weakly_connected_components(G):
        sub = G.subgraph(c)
        sub_order = [n for n in order if n in c]
        anc = {}
        for n in sub_order:
            s = set()
            for p in sub.predecessors(n):
                s.add(p); s |= anc[p]
            anc[n] = s
        desc = {}
        for n in reversed(sub_order):
            s = set()
            for q in sub.successors(n):
                s.add(q); s |= desc[q]
            desc[n] = s
        for n in c:
            n_older[n] = len(anc[n]); n_younger[n] = len(desc[n])
    nodes = sorted(G.nodes)
    return pd.DataFrame({
        "node": nodes,
        "depth": [depth[n] for n in nodes],
        "height": [height[n] for n in nodes],
        "n_older": [n_older[n] for n in nodes],
        "n_younger": [n_younger[n] for n in nodes],
        "component": [comp[n] for n in nodes],
        "component_size": [csize[n] for n in nodes],
    }).set_index("node")


def hasse_edges(G: nx.DiGraph) -> nx.DiGraph:
    """Transitive reduction: the minimal edge set implying the same ordering."""
    R = nx.transitive_reduction(G)
    R.add_edges_from((u, v, G.edges[u, v]) for u, v in R.edges)
    return R
