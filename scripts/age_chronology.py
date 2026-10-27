"""Formation crater times and model ages for craters >= 10 km, with and without
superposition constraints.

  python scripts/age_chronology.py --out outputs [--dmin 3] [--eps 0] [--tag main]

Reads outputs/interior_counts.parquet, nodes.csv.gz and edges.csv; writes
outputs/crater_ages_<tag>.parquet and prints a summary.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.ages import t_of_n1, epoch_boundaries_tau
from mcgraph.chronology import tau_grid, log_lik, fit_prior, gibbs, gibbs_arrival

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--dmin", default="3")
ap.add_argument("--eps", type=float, default=0.0)
ap.add_argument("--tag", default="main")
ap.add_argument("--sweeps", type=int, default=4000)
ap.add_argument("--prior", choices=["global", "bysize"], default="global")
ap.add_argument("--model", choices=["arrival", "independent"], default="arrival")
a = ap.parse_args()

r = pd.read_parquet(os.path.join(a.out, "interior_counts.parquet"))
e = pd.read_csv(os.path.join(a.out, "edges.csv"))
e = e[e.kept]
hosts = r.index.to_numpy()
pos = {h: i for i, h in enumerate(hosts)}

# order relations among hosts, through any path in the full graph
G = nx.DiGraph()
G.add_edges_from(zip(e.older, e.younger))
rel = set()
for h in hosts:
    if h in G:
        for d in nx.descendants(G, h):
            if d in pos:
                rel.add((h, d))
HG = nx.DiGraph()
HG.add_edges_from(rel)
HR = nx.transitive_reduction(HG)
cons = np.array([(pos[o], pos[y]) for o, y in HR.edges()])
print("host relations", len(rel), "after reduction", len(cons))

grid = tau_grid()
k = r[f"k_{a.dmin}"].to_numpy()
L = r[f"L_{a.dmin}"].to_numpy()
LL = log_lik(k, L, grid, eps=a.eps)

# empirical-Bayes prior by host size class
bins = [10, 20, 40, 80, 1e9] if a.prior == "bysize" else [0, 1e9]
cls = np.digitize(r.D_km.to_numpy(), bins) - 1
logprior = np.zeros_like(LL)
priors = {}
for c in range(len(bins) - 1):
    m = cls == c
    g = fit_prior(LL[m])
    priors[f"{bins[c]}-{bins[c+1]}"] = g.tolist()
    logprior[m] = np.log(g + 1e-300)[None, :]
lp = LL + logprior


def q_from(lp_rows, qs=(0.16, 0.5, 0.84)):
    P = np.exp(lp_rows - lp_rows.max(axis=1, keepdims=True))
    C = np.cumsum(P, axis=1); C /= C[:, -1:]
    return [grid[np.clip((C < q).sum(axis=1), 0, len(grid) - 1)] for q in qs]


u16, u50, u84 = q_from(lp)
c16, c50, c84 = u16.copy(), u50.copy(), u84.copy()
n_old = np.zeros(len(hosts), int); n_young = np.zeros(len(hosts), int)
for o, y in cons:
    n_young[o] += 1; n_old[y] += 1

# Gibbs on the constrained craters only
sub = np.unique(cons.ravel())
smap = {v: i for i, v in enumerate(sub)}
scons = np.array([(smap[o], smap[y]) for o, y in cons])
t0 = time.time()
if a.model == "arrival":
    S = gibbs_arrival(LL[sub], logprior[sub], scons, grid, n_sweeps=a.sweeps, burn=600, seed=1)
else:
    S = gibbs(lp[sub], scons, n_sweeps=a.sweeps, burn=600, seed=1)
np.savez_compressed(os.path.join(a.out, f"crater_ages_{a.tag}_samples.npz"),
                    nodes=hosts[sub], idx=S[::4], grid=grid)
print("gibbs", len(sub), "craters", round(time.time() - t0), "s")
ts = grid[S]
c16[sub], c50[sub], c84[sub] = np.quantile(ts, [0.16, 0.5, 0.84], axis=0)
# probability that each constrained crater's interior count is "too young" for its
# constrained age: posterior predictive P(k' <= k)
lam = ts * L[sub][None, :]
from scipy.stats import poisson
ppp = poisson.cdf(k[sub][None, :], lam).mean(axis=0)
p_low = np.full(len(hosts), np.nan); p_low[sub] = ppp

out = pd.DataFrame({
    "CRATER_ID": r.CRATER_ID.to_numpy(), "D_km": r.D_km.to_numpy(),
    "k": k, "Lambda": L, "n_older": n_old, "n_younger": n_young,
    "tau16_u": u16, "tau50_u": u50, "tau84_u": u84,
    "tau16": c16, "tau50": c50, "tau84": c84,
    "T16_u": t_of_n1(u16), "T50_u": t_of_n1(u50), "T84_u": t_of_n1(u84),
    "T16": t_of_n1(c16), "T50": t_of_n1(c50), "T84": t_of_n1(c84),
    "p_count_low": p_low,
}, index=pd.Index(hosts, name="node"))
out.to_parquet(os.path.join(a.out, f"crater_ages_{a.tag}.parquet"))
json.dump({"grid": grid.tolist(), "priors": priors, "dmin": a.dmin, "eps": a.eps, "model": a.model},
          open(os.path.join(a.out, f"crater_ages_{a.tag}_prior.json"), "w"))

con = out.n_older + out.n_younger > 0
w_u = np.log(out.tau84_u / out.tau16_u); w_c = np.log(out.tau84 / out.tau16)
wt_u = out.T84_u - out.T16_u; wt_c = out.T84 - out.T16
print(f"constrained craters: {con.sum()} of {len(out)}")
print("median width ratio (log tau):", round(float(np.median((w_c / w_u)[con])), 3))
print("median 68% width, Ga: unconstrained", round(float(wt_u[con].median()), 3),
      "constrained", round(float(wt_c[con].median()), 3))
for f in (0.8, 0.5):
    print(f"constrained craters with tau interval narrowed to <= {f:.0%}:",
          int(((w_c / w_u) <= f)[con].sum()))
shift = np.log(out.tau50 / out.tau50_u)
print("median |shift| in log tau:", round(float(shift[con].abs().median()), 3),
      "; >= x1.5:", int((shift[con].abs() >= np.log(1.5)).sum()))
print("interior counts low for constrained age (p<0.05):", int((out.p_count_low < 0.05).sum()))
# information: 68% width relative to the prior's
g = np.array(next(iter(priors.values()))); cg = np.cumsum(g)
pw = np.log(grid[np.searchsorted(cg, 0.84)] / grid[np.searchsorted(cg, 0.16)])
for f in (0.5, 0.33):
    iu = (w_u <= f * pw); ic = (w_c <= f * pw)
    print(f"craters with 68% interval <= {f:.0%} of the prior's: counts only {int(iu.sum())}, with order {int(ic.sum())}")
