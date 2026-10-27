"""Compare simulated superposition graphs with the observed graph on Noachian units.

  python scripts/sim_compare.py --out outputs --sims outputs/sim/seed1 outputs/sim/seed2 ...

The simulated patch has the area of the Noachian units of Tanaka et al. (2014). For
each simulation the snapshot is chosen where surviving craters >= 16 km match the
observed Noachian density. Three record filters are then applied before the graph
is built, standing in for erosion, infill, and catalog incompleteness:
  none          every crater that survived the impacts
  thin_random   each size bin thinned at random to the observed Noachian density
  thin_oldest   each size bin thinned to the observed density by removing the oldest
                craters first
and three edge sets: containment only; containment plus every rim overlap; and
containment plus 59% of rim overlaps (the fraction the pipeline decides).
Writes outputs/sim_compare_<erasure>.json.
"""
import argparse, json, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--sims", nargs="+", required=True)
ap.add_argument("--erasure", choices=["rim", "ejecta", "arc"], default="ejecta")
a = ap.parse_args()

AREA_N = 6.827976e7
BINS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 362.0001]

# ---------------- observed, Noachian units ----------------
nodes = pd.read_csv(os.path.join(a.out, "nodes.csv.gz"), low_memory=False)
edges = pd.read_csv(os.path.join(a.out, "edges.csv"))
edges = edges[edges.kept]
NN = len(nodes)
isN = (nodes.unit_period == "N").to_numpy()


def depth_from(o, y, n):
    d = np.zeros(n, dtype=np.int64)
    for _ in range(64):
        new = d.copy()
        np.maximum.at(new, y, d[o] + 1)
        if np.array_equal(new, d):
            break
        d = new
    return d


def summarize(o, y, n, mask):
    d = depth_from(o, y, n)
    indeg = np.bincount(y, minlength=n); outdeg = np.bincount(o, minlength=n)
    linked = mask & ((indeg + outdeg) > 0)
    dc = np.bincount(d[linked], minlength=8)[:8]
    ic = np.bincount(indeg[linked], minlength=8)[:8]
    ge = [int(dc[k:].sum()) for k in range(8)]
    return {"n_record": int(mask.sum()), "n_linked": int(linked.sum()),
            "depth_counts": dc.tolist(), "indeg_counts": ic.tolist(),
            "ratio_ge2_ge1": ge[2] / max(ge[1], 1), "ratio_ge3_ge2": ge[3] / max(ge[2], 1),
            "ratio_ge4_ge3": ge[4] / max(ge[3], 1), "max_depth": int(np.nonzero(dc)[0].max())}


obs_D = nodes.D_km.to_numpy()[isN]
obs_bin_density = np.histogram(obs_D, BINS)[0] / AREA_N
Dn = nodes.D_km.to_numpy()
edges = edges[(Dn[edges.older] <= 362) & (Dn[edges.younger] <= 362)]  # simulated size range
isN = isN & (Dn <= 362)
obs_D = Dn[isN]
obs_bin_density = np.histogram(obs_D, BINS)[0] / AREA_N
cont = (edges.method == "containment").to_numpy()
o_all, y_all = edges.older.to_numpy(), edges.younger.to_numpy()
observed = {
    "containment": summarize(o_all[cont], y_all[cont], NN, isN),
    "contain+decided_overlap": summarize(o_all, y_all, NN, isN),
    "N>=16": float((obs_D >= 16).sum() / AREA_N),
    "N>=1": float((obs_D >= 1).sum() / AREA_N),
    "coverage_<=362km": float((np.pi / 4 * obs_D[obs_D <= 362] ** 2).sum() / AREA_N),
    "bin_density": obs_bin_density.tolist(),
}
print("observed", {k: (v["depth_counts"][:6], round(v["ratio_ge2_ge1"], 3)) for k, v in observed.items()
                   if isinstance(v, dict)})

# ---------------- simulations ----------------
results = []
for sdir in a.sims:
    z = np.load(os.path.join(sdir, "sim_state.npz"))
    D, older, younger, contain = z["D"], z["older"], z["younger"], z["contain"]
    er = {"rim": z["erase_rim"], "ejecta": z["erase_ej"], "arc": z["erase_arc"]}[a.erasure]
    L = float(z["side"]); area = L * L
    n = len(D)
    idx = np.arange(n)
    big16 = D >= 16
    # snapshot matching N(>=16)
    grid = np.linspace(0.05, 1, 400) * n
    n16 = np.array([((big16[:int(g)]) & (er[:int(g)] >= int(g))).sum() / area for g in grid])
    ns = int(grid[np.argmin(np.abs(n16 - observed["N>=16"]))])
    alive0 = (idx < ns) & (er >= ns)
    rng = np.random.default_rng(7)
    bin_of = np.digitize(D, BINS) - 1
    filters = {"none": alive0.copy(), "thin_random": alive0.copy(), "thin_oldest": alive0.copy()}
    for b in range(len(BINS) - 1):
        members = np.nonzero(alive0 & (bin_of == b))[0]
        target = int(round(obs_bin_density[b] * area))
        if len(members) > target:
            drop = rng.choice(members, len(members) - target, replace=False)
            filters["thin_random"][drop] = False
            filters["thin_oldest"][np.sort(members)[: len(members) - target]] = False  # oldest = lowest index
    dec = np.random.default_rng(11).random(len(older)) < 0.59
    # decided fraction of overlap pairs by the smaller crater's diameter, as observed
    ds = np.minimum(D[older], D[younger])
    p_dec = np.interp(np.log(ds), np.log([1.4, 2.8, 5.7, 11.3, 32, 160]),
                      [0.69, 0.59, 0.52, 0.44, 0.46, 0.39])
    dec_size = np.random.default_rng(12).random(len(older)) < p_dec
    for fname, alive in filters.items():
        Da = D[alive]
        for ename, emask in (("containment", contain), ("contain+all_overlap", np.ones_like(contain)),
                             ("contain+59%_overlap", contain | dec),
                             ("contain+sized_overlap", contain | dec_size)):
            m = alive[older] & alive[younger] & emask
            s = summarize(older[m], younger[m], n, alive)
            s |= {"sim": os.path.basename(sdir), "tau": ns / area, "filter": fname, "edges": ename,
                  "N>=1": float((Da >= 1).sum() / area), "N>=16": float((Da >= 16).sum() / area),
                  "coverage": float((np.pi / 4 * Da ** 2).sum() / area)}
            results.append(s)
            print(f"{s['sim']} tau={s['tau']:.4f} {fname:11s} {ename:20s} depth {s['depth_counts'][:7]} "
                  f"ge2/ge1 {s['ratio_ge2_ge1']:.3f} cov {s['coverage']:.3f}")
json.dump({"observed_noachian": observed, "simulated": results, "erasure": a.erasure},
          open(os.path.join(a.out, f"sim_compare_{a.erasure}.json"), "w"), indent=1)
