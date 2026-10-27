"""Monte Carlo cratering of a flat periodic patch, and the superposition graph its
surviving craters would produce.

  python scripts/sim_cratering.py --out outputs/sim --side 4000 --tau-max 0.04 --seed 1

Craters >= 1 km arrive one at a time at uniform random positions, with diameters
drawn from the Ivanov (2001) production function (1-362 km, its validity range).
Crater time tau is the number of craters >= 1 km per km^2 produced so far.

Erasure (a crater is removed from the record and never returns):
  rim      a later, larger crater whose rim encloses the older crater's center
  ejecta   (optional) a later crater's continuous ejecta, thickness
           t = 0.14 R^0.74 (r/R)^-3 m (McGetchin et al., 1973), at least as thick as
           the older crater is deep; depth = min(0.2 D, 0.36 D^0.49) km
  arc      (optional) later craters, of any size, have together cut away at least
           --arc-frac of the older crater's rim

The graph uses the pipeline's geometric rules on the craters that survive at each
snapshot: containment (d + r_small <= r_large) and rim overlap (penetration
>= 0.1 r_small), with edges from older to younger. Ejecta-zone and touching pairs
are not linked. Results are written per snapshot and per scenario.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.ages import pf_rel

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--side", type=float, default=4000.0, help="patch side, km")
ap.add_argument("--tau-max", type=float, default=0.04)
ap.add_argument("--dmax", type=float, default=362.0)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--arc-frac", type=float, default=0.5, help="rim fraction lost that removes a crater")
ap.add_argument("--snapshots", default="0.0025,0.005,0.01,0.015,0.02,0.03,0.04")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng(a.seed)
L = a.side
area = L * L
n = int(round(a.tau_max * area))

# diameters by inverse transform of N(>=D)/N(>=1)
dg = np.geomspace(1.0, a.dmax, 4000)
cg = 1.0 - pf_rel(dg) / pf_rel(1.0)          # CDF of D among craters >= 1 km
cg /= cg[-1]
D = np.interp(rng.random(n), cg, dg)
r = D / 2
xy = rng.random((n, 2)) * L
tau = (np.arange(n) + 1) / area
print(f"{n:,} craters on {L:.0f} km patch to tau = {a.tau_max}")

tree = cKDTree(xy, boxsize=L)
depth_km = np.minimum(0.2 * D, 0.36 * D ** 0.49)


def neighbors(radius, idx):
    """Pairs (j, i) with i within radius[j] of j, for j in idx (periodic)."""
    out_j, out_i = [], []
    for chunk in np.array_split(idx, max(1, len(idx) // 200000)):
        lists = tree.query_ball_point(xy[chunk], radius[chunk], workers=-1)
        lens = np.fromiter((len(l) for l in lists), int, len(lists))
        out_j.append(np.repeat(chunk, lens))
        out_i.append(np.fromiter((v for l in lists for v in l), int, lens.sum()))
    return np.concatenate(out_j), np.concatenate(out_i)


def pdist(i, j):
    d = np.abs(xy[i] - xy[j])
    d = np.minimum(d, L - d)
    return np.hypot(d[:, 0], d[:, 1])


t0 = time.time()
# erasure: for every crater j, older craters i (i < j) it removes
jj, ii = neighbors(3 * r, np.arange(n))
keep = ii < jj
jj, ii = jj[keep], ii[keep]
d = pdist(ii, jj)
by_rim = (d < r[jj]) & (r[ii] < r[jj])
R_m = r[jj] * 1000
with np.errstate(divide="ignore"):
    thick_km = 0.14 * R_m ** 0.74 * np.clip(d / r[jj], 1, None) ** -3 / 1000
by_ejecta = (d >= r[jj]) & (thick_km >= depth_km[ii])
erase_rim = np.full(n, n, dtype=np.int64)
np.minimum.at(erase_rim, ii[by_rim], jj[by_rim])
erase_ej = erase_rim.copy()
np.minimum.at(erase_ej, ii[by_ejecta], jj[by_ejecta])
print("erasure pairs", int(by_rim.sum()), int(by_ejecta.sum()), round(time.time() - t0), "s")

# candidate graph pairs: each crater with smaller craters within 2 r of it
gj, gi = neighbors(2 * r, np.arange(n))
m = (gi != gj) & ((r[gi] < r[gj]) | ((r[gi] == r[gj]) & (gi > gj)))
big, small = gj[m], gi[m]
dd = pdist(big, small)
contain = dd + r[small] <= r[big]
overlap = ~contain & ((r[big] + r[small] - dd) >= 0.1 * r[small])

# rim-arc erasure: an older crater is lost from the record once later craters have cut
# away a fraction >= arc_frac of its rim (64 sample points on the rim)
x_ = dd < r[big] + r[small]
po, py = np.minimum(big[x_], small[x_]), np.maximum(big[x_], small[x_])
order = np.lexsort((py, po)); po, py = po[order], py[order]
th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
bits = np.zeros(len(po), dtype=np.uint64)
for c in range(0, len(po), 500000):
    s_ = slice(c, c + 500000)
    dx = xy[py[s_]] - xy[po[s_]]
    dx -= L * np.round(dx / L)
    px = r[po[s_]][:, None] * np.cos(th)[None, :] - dx[:, 0:1]
    pyy = r[po[s_]][:, None] * np.sin(th)[None, :] - dx[:, 1:2]
    inside = (px ** 2 + pyy ** 2) < r[py[s_]][:, None] ** 2
    bits[s_] = (inside.astype(np.uint64) << np.arange(64, dtype=np.uint64)[None, :]).sum(axis=1, dtype=np.uint64)
k = 1
while k < len(po):
    same = po[k:] == po[:-k]
    if not same.any():
        break
    new = bits.copy()
    new[k:] |= np.where(same, bits[:-k], np.uint64(0))
    bits = new
    k *= 2
hit = np.bitwise_count(bits) >= int(round(a.arc_frac * 64))
erase_arc_only = np.full(n, n, dtype=np.int64)
np.minimum.at(erase_arc_only, po[hit], py[hit])
erase_arc = np.minimum(erase_ej, erase_arc_only)
print("rim-arc erasures", int((erase_arc < erase_ej).sum()), round(time.time() - t0), "s")
sel = contain | overlap
big, small, contain = big[sel], small[sel], contain[sel]
older = np.minimum(big, small); younger = np.maximum(big, small)
print("graph candidate pairs", len(big), round(time.time() - t0), "s")


def graph_stats(alive, decided_frac=1.0, seed=0):
    """Depth and in-degree of surviving craters given an alive mask."""
    e = alive[older] & alive[younger]
    if decided_frac < 1:
        rr = np.random.default_rng(seed).random(len(older))
        e &= contain | (rr < decided_frac)
    o, y = older[e], younger[e]
    depth = np.zeros(n, dtype=np.int64)
    for _ in range(64):  # longest-path depth by relaxation; converges in max-depth + 1 passes
        new = depth.copy()
        np.maximum.at(new, y, depth[o] + 1)
        if np.array_equal(new, depth):
            break
        depth = new
    indeg = np.bincount(y, minlength=n)
    outdeg = np.bincount(o, minlength=n)
    linked = alive & ((indeg + outdeg) > 0)
    dc = np.bincount(depth[linked], minlength=8)
    ic = np.bincount(indeg[linked], minlength=8)
    return {"n_alive": int(alive.sum()), "n_linked": int(linked.sum()), "n_edges": int(e.sum()),
            "depth_counts": dc.tolist(), "indeg_counts": ic.tolist(),
            "max_depth": int(depth[linked].max()) if linked.any() else 0}


def density(alive, ns):
    Da = D[:ns][alive[:ns]]
    return {f"N>={k}": float((Da >= k).sum() / area) for k in (1, 2, 4, 8, 16, 32, 64)} | {
        "coverage": float((np.pi / 4 * Da ** 2).sum() / area)}


obs_dens = None
res = []
for ts in map(float, a.snapshots.split(",")):
    ns = int(round(ts * area))
    for scen, er in (("rim", erase_rim), ("rim+ejecta", erase_ej), ("rim+ejecta+arc", erase_arc)):
        alive = np.zeros(n, bool)
        alive[:ns] = er[:ns] >= ns
        row = {"tau": ts, "scenario": scen, "filter": "none"} | density(alive, ns)
        row |= {f"all_{k}": v for k, v in graph_stats(alive).items()}
        row |= {f"dec59_{k}": v for k, v in graph_stats(alive, 0.59).items()}
        res.append(row)
        print(ts, scen, row["all_n_alive"], row["all_depth_counts"][:6],
              row["dec59_depth_counts"][:6], round(row["coverage"], 3))
ej = (d >= r[jj]) & (thick_km > 0.001)  # ejecta-zone pairs, for depth-dependent burial later
np.savez_compressed(os.path.join(a.out, "sim_state.npz"), D=D, xy=xy, erase_rim=erase_rim,
                    erase_ej=erase_ej, erase_arc=erase_arc, erase_arc_only=erase_arc_only,
                    ej_i=ii[ej].astype(np.int32), ej_j=jj[ej].astype(np.int32),
                    ej_thick=thick_km[ej].astype(np.float32), older=older, younger=younger, contain=contain, side=L)
pd.DataFrame(res).to_json(os.path.join(a.out, "sim_snapshots.json"), orient="records", indent=1)
print("done", round(time.time() - t0), "s")
