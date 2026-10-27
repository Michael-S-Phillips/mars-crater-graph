"""Do interior crater densities agree with superposition order?

  python scripts/age_pair_test.py --out outputs

For every graph edge between two craters >= 10 km, the counts on the two interiors
are compared with the conditional binomial: given k = k_old + k_young, k_old ~
Bin(k, rho L_old / (rho L_old + L_young)), where rho = tau_old / tau_young.

1. Per diameter threshold, log rho across edges is modeled as Normal(mu, sigma),
   with counting noise integrated out; the fraction of edges whose older crater
   has the sparser interior is Phi(-mu / sigma).
2. Small (1-2 km) and large (>= 2 km) craters are disjoint samples of each interior.
   A bivariate normal for (log rho_small, log rho_large) gives the fraction of
   edges whose order is inverted by small craters but not by large ones.
Bootstrap over edges gives 95% intervals. Writes outputs/age_pair_test.json.
"""
import argparse, json, os
import numpy as np
import pandas as pd
from scipy.stats import binom, norm, multivariate_normal
from scipy.optimize import minimize

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--boot", type=int, default=200)
a = ap.parse_args()
r = pd.read_parquet(os.path.join(a.out, "interior_counts.parquet"))
e = pd.read_csv(os.path.join(a.out, "edges.csv"))
e = e[e.kept]
E = e[e.older.isin(r.index) & e.younger.isin(r.index)].reset_index(drop=True)
lg = np.linspace(-5, 5, 201); dl = lg[1] - lg[0]


# the younger crater of an edge is often itself counted in the older crater's interior;
# it defines the edge, so it is removed from the older crater's count
mem = pd.read_parquet(os.path.join(a.out, "interior_members.parquet"))
inmem = set(zip(mem.host, mem.crater))
self_in = np.array([(o, y) in inmem for o, y in zip(E.older, E.younger)])
Dy = r.loc[E.younger, "D_km"].to_numpy()
print("edges whose younger crater is counted in the older interior:", int(self_in.sum()))


def arr(col):
    ko, ky = r.loc[E.older, col].to_numpy().copy(), r.loc[E.younger, col].to_numpy()
    if col.startswith("k_"):
        ko = ko - (self_in & (Dy >= float(col[2:]))).astype(int)
    return ko, ky


def lik_matrix(ko, ky, Lo, Ly):
    K = ko + ky
    q = (np.exp(lg)[None, :] * Lo[:, None]) / (np.exp(lg)[None, :] * Lo[:, None] + Ly[:, None])
    return binom.pmf(ko[:, None], K[:, None], q)


def fit1(M):
    def nll(p):
        w = norm.pdf(lg, p[0], np.exp(p[1])) * dl
        return -np.sum(np.log(M @ w + 1e-300))
    res = minimize(nll, [0.6, np.log(0.4)], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-4})
    return res.x[0], np.exp(res.x[1])


def usable(ko, ky, Lo, Ly):
    return (ko + ky > 0) & (Lo > 0) & (Ly > 0)


out = {"n_edges": int(len(E)), "n_self_counted": int(self_in.sum()), "by_dmin": {}, "by_method": {}}
rng = np.random.default_rng(0)
for dm in ["1", "1.5", "2", "3"]:
    (ko, ky), (Lo, Ly) = arr(f"k_{dm}"), arr(f"L_{dm}")
    u = usable(ko, ky, Lo, Ly)
    M = lik_matrix(ko[u], ky[u], Lo[u], Ly[u])
    mu, sd = fit1(M)
    bs = []
    for _ in range(a.boot):
        i = rng.integers(0, M.shape[0], M.shape[0])
        m_, s_ = fit1(M[i]); bs.append((np.exp(m_), norm.cdf(-m_ / s_)))
    bs = np.array(bs)
    out["by_dmin"][dm] = {"n": int(u.sum()), "ratio": float(np.exp(mu)), "sigma": float(sd),
                          "ratio_95": np.quantile(bs[:, 0], [0.025, 0.975]).tolist(),
                          "frac_inverted": float(norm.cdf(-mu / sd)),
                          "frac_inverted_95": np.quantile(bs[:, 1], [0.025, 0.975]).tolist()}
    print(dm, out["by_dmin"][dm])
    for meth in ["containment", "arc", "topo+arc"]:
        s = u & (E.method.to_numpy() == meth)
        m_, s_ = fit1(lik_matrix(ko[s], ky[s], Lo[s], Ly[s]))
        out["by_method"].setdefault(meth, {})[dm] = {"n": int(s.sum()), "ratio": float(np.exp(m_)),
                                                    "frac_inverted": float(norm.cdf(-m_ / s_))}

# disjoint small (1-2 km) and large (>= 2 km) samples
(k1o, k1y), (k2o, k2y) = arr("k_1"), arr("k_2")
(L1o, L1y), (L2o, L2y) = arr("L_1"), arr("L_2")
kso, ksy, Lso, Lsy = k1o - k2o, k1y - k2y, L1o - L2o, L1y - L2y
u = usable(kso, ksy, Lso, Lsy) & usable(k2o, k2y, L2o, L2y)
g2 = np.linspace(-4, 4, 101); d2 = g2[1] - g2[0]


def lm(ko, ky, Lo, Ly):
    K = ko + ky
    q = (np.exp(g2)[None, :] * Lo[:, None]) / (np.exp(g2)[None, :] * Lo[:, None] + Ly[:, None])
    return binom.pmf(ko[:, None], K[:, None], q)


MS, MLg = lm(kso[u], ksy[u], Lso[u], Lsy[u]), lm(k2o[u], k2y[u], L2o[u], L2y[u])
X, Y = np.meshgrid(g2, g2, indexing="ij")
pts = np.stack([X.ravel(), Y.ravel()], 1)


def dens(p):
    m1, m2, l1, l2, z = p
    s1, s2, rho = np.exp(l1), np.exp(l2), np.tanh(z)
    cov = [[s1 * s1, rho * s1 * s2], [rho * s1 * s2, s2 * s2]]
    return multivariate_normal([m1, m2], cov).pdf(pts).reshape(X.shape) * d2 * d2


def nll2(p, MS=MS, MLg=MLg):
    W = dens(p)
    return -np.sum(np.log(np.einsum("ni,ij,nj->n", MS, W, MLg) + 1e-300))


res = minimize(nll2, [0.5, 0.7, np.log(0.4), np.log(0.4), 0.5], method="Nelder-Mead",
               options={"maxiter": 4000, "xatol": 1e-4, "fatol": 1e-4})
W = dens(res.x); W /= W.sum()
f_small_inv = float(W[X < 0].sum()); f_large_inv = float(W[Y < 0].sum())
f_small_only = float(W[(X < 0) & (Y > 0)].sum())
out["small_vs_large"] = {"n": int(u.sum()), "ratio_small": float(np.exp(res.x[0])),
                         "ratio_large": float(np.exp(res.x[1])),
                         "sigma_small": float(np.exp(res.x[2])), "sigma_large": float(np.exp(res.x[3])),
                         "corr": float(np.tanh(res.x[4])),
                         "frac_inverted_small": f_small_inv, "frac_inverted_large": f_large_inv,
                         "frac_inverted_small_only": f_small_only}
print(out["small_vs_large"])
json.dump(out, open(os.path.join(a.out, "age_pair_test.json"), "w"), indent=1)
