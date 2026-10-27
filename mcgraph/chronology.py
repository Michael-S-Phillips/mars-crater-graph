"""Bayesian crater chronology constrained by superposition order.

Each host crater i has a formation crater time tau_i (N(>=1 km) per km^2 since it
formed; see ages.py). Its own interior count gives a likelihood,

    k_i ~ Poisson(tau_i * Lambda_i),

optionally mixed with a "resurfaced interior" term: with probability eps the
interior was reset at a uniformly random crater time after the crater formed,
so that the count records only part of tau_i:

    p(k | tau) = (1 - eps) Pois(k; tau Lambda) + eps * P(k + 1, tau Lambda) / (tau Lambda),

where P is the regularized lower incomplete gamma function.

The prior for tau_i is the population distribution of formation crater times for
craters of similar size, estimated from all host craters by empirical Bayes
(EM on a log-spaced grid). Superposition edges add hard constraints
tau_older >= tau_younger. Because the constraints only couple craters within a
connected cluster of the graph, unconstrained craters keep their closed-form
posterior, and the constrained ones are sampled by Gibbs sampling on the grid,
updating one color class of the constraint graph at a time.

This mirrors Bayesian chronological modeling of stratigraphic sequences
(Buck et al., 1991; Bronk Ramsey, 2009) with crater counts in place of
radiocarbon determinations.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln, gammainc


def tau_grid(n=400, lo=1e-5, hi=0.15):
    return np.geomspace(lo, hi, n)


def log_lik(k, L, grid, eps=0.0):
    """Log-likelihood of each crater's count on the tau grid (rows: craters)."""
    k = np.asarray(k, float)[:, None]
    mu = grid[None, :] * np.asarray(L, float)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        lp = np.where(mu > 0, k * np.log(np.where(mu > 0, mu, 1.0)) - mu - gammaln(k + 1),
                      np.where(k == 0, 0.0, -np.inf))
        if eps <= 0:
            return lp
        res = np.where(mu > 1e-12, gammainc(k + 1, mu) / np.where(mu > 1e-12, mu, 1.0),
                       np.where(k == 0, 1.0, 0.0))
        return np.log((1 - eps) * np.exp(lp) + eps * res + 1e-300)


def fit_prior(LL, iters=400, smooth=2.0):
    """Empirical-Bayes population distribution on the grid (EM, lightly smoothed)."""
    G = LL.shape[1]
    g = np.full(G, 1.0 / G)
    Lmax = LL.max(axis=1, keepdims=True)
    W = np.exp(LL - Lmax)
    x = np.arange(G)
    kern = np.exp(-0.5 * ((x[:, None] - x[None, :]) / smooth) ** 2) if smooth else None
    for _ in range(iters):
        R = W * g[None, :]
        R /= R.sum(axis=1, keepdims=True)
        g = R.mean(axis=0)
        if kern is not None:
            g = kern @ g
            g /= g.sum()
    return g


def _bounds_index(C, lo, hi, u):
    """Sample grid indices in [lo, hi] from row-wise CDFs C by inversion."""
    n = len(lo)
    rows = np.arange(n)
    c_lo = np.where(lo > 0, C[rows, np.maximum(lo - 1, 0)], 0.0)
    c_hi = C[rows, hi]
    target = c_lo + u * (c_hi - c_lo)
    idx = (C < target[:, None]).sum(axis=1)
    return np.clip(idx, lo, hi)


def gibbs(logpost, edges, n_sweeps=3000, burn=500, seed=0, colors=None):
    """Gibbs sampler for grid indices under order constraints.

    logpost  (n, G) unnormalized log posterior of each node without constraints
    edges    (m, 2) array of (older, younger) node indices (0..n-1)
    Returns  samples (n_sweeps, n) of grid indices.
    """
    rng = np.random.default_rng(seed)
    n, G = logpost.shape
    P = np.exp(logpost - logpost.max(axis=1, keepdims=True))
    C = np.cumsum(P, axis=1)
    C /= C[:, -1:]
    younger_of = [[] for _ in range(n)]
    older_of = [[] for _ in range(n)]
    for o, y in edges:
        younger_of[o].append(y)
        older_of[y].append(o)
    # padded neighbor arrays for vectorized min/max
    def pad(lists, fill):
        w = max(1, max(len(l) for l in lists))
        a = np.full((n, w), fill, dtype=np.int64)
        for i, l in enumerate(lists):
            a[i, :len(l)] = l
        return a
    Y = pad(younger_of, n)   # index n -> sentinel
    O = pad(older_of, n)
    if colors is None:
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(n))
        g.add_edges_from(map(tuple, edges))
        col = nx.greedy_color(g, strategy="largest_first")
        colors = np.array([col[i] for i in range(n)])
    classes = [np.where(colors == c)[0] for c in np.unique(colors)]
    # feasible start: posterior medians, then push older nodes up in topological order
    idx = (C < 0.5).sum(axis=1)
    import networkx as nx
    dg = nx.DiGraph()
    dg.add_nodes_from(range(n))
    dg.add_edges_from(map(tuple, edges))
    for v in reversed(list(nx.topological_sort(dg))):  # youngest first
        ys = younger_of[v]
        if ys:
            idx[v] = max(idx[v], max(idx[y] for y in ys))
    out = np.empty((n_sweeps, n), dtype=np.int16)
    ext = np.concatenate([idx, [0]])  # sentinel slot for padding
    for s in range(burn + n_sweeps):
        for cl in classes:
            ext_lo = ext.copy(); ext_lo[n] = 0
            lo = ext_lo[Y[cl]].max(axis=1)
            ext_hi = ext.copy(); ext_hi[n] = G - 1
            hi = ext_hi[O[cl]].min(axis=1)
            hi = np.maximum(hi, lo)
            ext[cl] = _bounds_index(C[cl], lo, hi, rng.random(len(cl)))
        if s >= burn:
            out[s - burn] = ext[:n]
    return out


def gibbs_arrival(LL, logg, edges, grid, n_sweeps=3000, burn=500, seed=0):
    """Gibbs sampler with the impact-arrival prior.

    Craters with no older neighbor (roots) take the population prior g. A crater
    that superposes others formed after all of them, and with impacts arriving at
    a constant rate in crater time it is uniform on [0, min tau_parents]. Its prior
    density 1/min(tau_parents) enters the parents' conditionals, which cancels the
    spurious push toward old ages that independent priors on many younger craters
    would cause (the same problem the uniform-phase prior solves in radiocarbon
    sequence models).

    LL     (n, G) log-likelihood on the grid
    logg   (G,) log population prior for roots
    edges  (m, 2) (older, younger) pairs, transitively reduced
    """
    rng = np.random.default_rng(seed)
    n, G = LL.shape
    parents = [[] for _ in range(n)]
    children = [[] for _ in range(n)]
    for o, y in edges:
        children[o].append(y)
        parents[y].append(o)
    root = np.array([len(p) == 0 for p in parents])
    n_single = np.array([sum(len(parents[c]) == 1 for c in ch) for ch in children])
    return _gibbs_arrival_core(LL, logg, parents, children, root, n_single, n_sweeps, burn, rng, grid)


def _gibbs_arrival_core(LL, logg, parents, children, root, n_single, n_sweeps, burn, rng, grid):
    import networkx as nx
    n, G = LL.shape
    lg = np.log(grid)
    # grid cells are log-spaced, so a density uniform in tau has cell mass proportional to the cell width
    ldt = np.log(np.gradient(grid))
    lgg = logg if np.ndim(logg) == 2 else np.broadcast_to(logg, LL.shape)
    base = (LL + np.where(root[:, None], lgg, ldt[None, :])
            - n_single[:, None] * lg[None, :])
    # (parent, child) pairs whose child has several parents: handled dynamically
    mp = [(p, c) for c in range(n) if len(parents[c]) > 1 for p in parents[c]]
    dyn = np.zeros(n, bool)
    for p, c in mp:
        dyn[p] = True
    # static CDFs
    P = np.exp(base - base.max(axis=1, keepdims=True))
    C = np.cumsum(P, axis=1); C /= C[:, -1:]
    def pad(lists, fill):
        w = max(1, max(len(l) for l in lists))
        a = np.full((n, w), fill, dtype=np.int64)
        for i, l in enumerate(lists):
            a[i, :len(l)] = l
        return a
    Y = pad(children, n)
    O = pad(parents, n)
    g = nx.Graph(); g.add_nodes_from(range(n))
    g.add_edges_from((o, y) for y in range(n) for o in parents[y])
    # color so that a node, its parents, its children and its co-parents differ
    for p, c in mp:
        for q in parents[c]:
            if q != p:
                g.add_edge(p, q)
    col = nx.greedy_color(g, strategy="largest_first")
    colors = np.array([col[i] for i in range(n)])
    classes = [np.where(colors == c)[0] for c in np.unique(colors)]
    mp_p = np.array([p for p, c in mp], dtype=np.int64)
    mp_c = np.array([c for p, c in mp], dtype=np.int64)
    # feasible start
    idx = (C < 0.5).sum(axis=1)
    dg = nx.DiGraph(); dg.add_nodes_from(range(n))
    dg.add_edges_from((o, y) for y in range(n) for o in parents[y])
    for v in reversed(list(nx.topological_sort(dg))):
        if children[v]:
            idx[v] = max(idx[v], max(idx[y] for y in children[v]))
    ext = np.concatenate([idx, [0]])
    out = np.empty((n_sweeps, n), dtype=np.int16)
    for s in range(burn + n_sweeps):
        for cl in classes:
            e_lo = ext.copy(); e_lo[n] = 0
            lo = e_lo[Y[cl]].max(axis=1)
            e_hi = ext.copy(); e_hi[n] = G - 1
            hi = np.maximum(e_hi[O[cl]].min(axis=1), lo)
            new = _bounds_index(C[cl], lo, hi, rng.random(len(cl)))
            d = dyn[cl]
            if d.any():
                dn = cl[d]
                loc = {v: i for i, v in enumerate(dn)}
                corr = np.zeros((len(dn), G))
                sel = np.isin(mp_p, dn)
                for p, c in zip(mp_p[sel], mp_c[sel]):
                    others = [q for q in parents[c] if q != p]
                    m = min(ext[q] for q in others)
                    # -log min(tau, tau_m): -log tau below the co-parent minimum, constant above
                    corr[loc[p]] += np.where(np.arange(G) <= m, -lg, -lg[m])
                bd = base[dn] + corr
                Pd = np.exp(bd - bd.max(axis=1, keepdims=True))
                Cd = np.cumsum(Pd, axis=1); Cd /= Cd[:, -1:]
                new[d] = _bounds_index(Cd, lo[d], hi[d], rng.random(len(dn)))
            ext[cl] = new
        if s >= burn:
            out[s - burn] = ext[:n]
    return out


def quantiles_from_logpost(logpost, grid, qs=(0.16, 0.5, 0.84)):
    P = np.exp(logpost - logpost.max(axis=1, keepdims=True))
    C = np.cumsum(P, axis=1)
    C /= C[:, -1:]
    return np.stack([grid[np.clip((C < q).sum(axis=1), 0, len(grid) - 1)] for q in qs], axis=1)
