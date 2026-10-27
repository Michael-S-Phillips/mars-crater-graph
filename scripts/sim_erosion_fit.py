"""Fit crater degradation rates to Noachian and Hesperian terrain with the cratering
simulation, then test the fit against the superposition graph.

  python scripts/sim_erosion_fit.py --out outputs --sims outputs/sim/seed1 ... [--tb 3.7] [--qmin 0.05]

Degradation model. A fresh crater has depth d0 = min(0.2 D, 0.28 D^0.51) km; the
second term is a power law through the 90th percentile of MOLA rim-to-floor depths of
rim-degradation-class 4 craters (Robbins and Hynek, 2012 depths). Depth then decreases
by erosion and infill at eps_early (m/Myr) before T_b and eps_late after it. A crater
leaves the record when its depth falls below qmin * d0. Ejecta from later impacts
buries an older crater when the ejecta is at least as thick as the older crater is
deep at that time, so degraded craters are buried more easily. Rim and rim-arc
erasure by later impacts are as in sim_cratering.py. Model ages use the Hartmann and
Neukum (2001) chronology; the surface's total crater time is set so that surviving
craters >= 16 km match the observed density.

Hesperian units: the same Noachian substrate is resurfaced at T_r by a deposit of
thickness h, which removes craters shallower than h and fills the rest by h;
cratering and degradation continue afterward.

Fitted observations: the size-frequency distribution in factor-of-2 bins (1-362 km)
and the distribution of relative depth q = depth/d0 for craters >= 5 km. The graph
statistics (chain depth, older-neighbor counts) are not fitted; they are predicted.
Writes outputs/sim_erosion_fit.json.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.ages import t_of_n1

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--sims", nargs="+", required=True)
ap.add_argument("--tb", type=float, default=3.7, help="Ga, end of the early degradation rate")
ap.add_argument("--qmin", type=float, default=0.05)
ap.add_argument("--tag", default="main")
ap.add_argument("--grid-sims", type=int, default=1, help="simulations used for the grid search")
ap.add_argument("--fit-min-km", type=float, default=2.0, help="smallest size bin in the fit")
ap.add_argument("--model", choices=["linear", "diffusion"], default="linear")
ap.add_argument("--thin12", type=int, default=1, help="1: thin 1-2 km craters at random to the observed count; 2: oldest first; 0: none")
ap.add_argument("--decision", choices=["size", "deg"], default="deg",
                help="overlap decision rates by the smaller crater's size, or by both craters' degradation class")
ap.add_argument("--best", default=None, help="eps_early,eps_late: skip the grid search")
ap.add_argument("--ee", default=None, help="lo,hi,n: geometric grid of early rates")
ap.add_argument("--el", default=None, help="lo,hi,n: geometric grid of late rates")
a = ap.parse_args()

AREA = {"N": 6.827976e7, "H": 6.195752e7}
BINS = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 362.0001])
QB = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 10])
TB = a.tb * 1000.0  # Myr


def d0_of(D):
    return np.minimum(0.2 * D, 0.28 * D ** 0.51)


# ---------------- observations ----------------
nodes = pd.read_csv(os.path.join(a.out, "nodes.csv.gz"), low_memory=False)
deps = pd.read_parquet(os.path.join(a.out, "crater_depths_matched.parquet"))
deps["q"] = deps.depth_km / d0_of(deps.D_km)
obs = {}
for u in ("N", "H"):
    D = nodes.D_km[(nodes.unit_period == u) & (nodes.D_km <= 362)].to_numpy()
    q = deps.q[(deps.unit_period == u) & (deps.D_km >= 5) & (deps.D_km <= 362)].to_numpy()
    obs[u] = {"sfd": np.histogram(D, BINS)[0].astype(float), "N16": float((D >= 16).sum()),
              "qhist": np.histogram(np.clip(q, 0, 9.99), QB)[0].astype(float)}


# observation model for rim-overlap decisions, from the real pairs: P(decided) by the
# rim degradation classes of the larger and smaller crater (smaller crater >= 3 km),
# and P(class | relative depth) to assign simulated craters a class
pairs = pd.read_parquet(os.path.join(a.out, "pairs_evidence.parquet"), columns=["big", "small", "kind", "confidence"])
ov = pairs[pairs.kind == "overlap"]
degn = nodes.DEG_RIM.to_numpy()
ok = (nodes.D_km.to_numpy()[ov.small] >= 3) & ~np.isnan(degn[ov.big]) & ~np.isnan(degn[ov.small])
dec_tab = (pd.DataFrame({"b": degn[ov.big][ok], "s": degn[ov.small][ok],
                         "d": ov.confidence.isin(["high", "medium"]).to_numpy()[ok]})
           .pivot_table(index="b", columns="s", values="d", aggfunc="mean").reindex(index=[1, 2, 3, 4], columns=[1, 2, 3, 4]).to_numpy())
QC = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 10])
dd3 = deps[(deps.D_km >= 3) & deps.DEG_RIM.notna()]
p_deg_q = np.array([np.bincount(dd3.DEG_RIM[(dd3.q >= QC[k]) & (dd3.q < QC[k + 1])].astype(int) - 1, minlength=4)
                    for k in range(len(QC) - 1)], float) + 0.5
p_deg_q /= p_deg_q.sum(1, keepdims=True)


def poisson_dev(o, e):
    e = np.maximum(e, 1e-9)
    t = np.where(o > 0, o * np.log(np.where(o > 0, o, 1) / e), 0.0)
    return 2 * float(np.sum(t - (o - e)))


def multinom_dev(o, p):
    p = np.maximum(p / p.sum(), 1e-9)
    e = o.sum() * p
    return 2 * float(np.sum(np.where(o > 0, o * np.log(np.where(o > 0, o, 1) / e), 0.0)))


# ---------------- simulations ----------------
sims = []
for sdir in a.sims:
    z = np.load(os.path.join(sdir, "sim_state.npz"))
    D = z["D"]; n = len(D); L = float(z["side"])
    # fresh depths scatter about 0.85 of the envelope (lognormal, sigma 0.15), so the
    # envelope sits near the 90th percentile of fresh depths, as it does for the data
    scatter = np.exp(np.log(0.85) + 0.15 * np.random.default_rng(3).standard_normal(n))
    s = {"name": os.path.basename(sdir), "D": D, "d0": d0_of(D) * scatter, "d0env": d0_of(D), "n": n, "area": L * L,
         "imp": np.minimum(z["erase_rim"], z["erase_arc_only"]),
         "ej_i": z["ej_i"].astype(np.int64), "ej_j": z["ej_j"].astype(np.int64), "ej_t": z["ej_thick"].astype(float),
         "older": z["older"], "younger": z["younger"], "contain": z["contain"],
         "u": np.random.default_rng(4).random(n)}
    ds = np.minimum(D[s["older"]], D[s["younger"]])
    p_dec = np.interp(np.log(ds), np.log([1.4, 2.8, 5.7, 11.3, 32, 160]), [0.69, 0.59, 0.52, 0.44, 0.46, 0.39])
    s["dec"] = s["contain"] | (np.random.default_rng(12).random(len(ds)) < p_dec)
    sims.append(s)
    print("loaded", s["name"], n)


def zcum(t_myr, ee, el):
    """Cumulative degradation between time t (Myr before present) and the present.
    linear: depth lost in km (rates in m/Myr); diffusion: kappa*t in m^2 (kappa in m^2/Myr)."""
    k = el * np.minimum(t_myr, TB) + ee * np.maximum(t_myr - TB, 0)
    return k / 1000.0 if a.model == "linear" else k


def depth_after(d0, D, dz):
    """Depth of a crater of fresh depth d0 (km), diameter D (km) after degradation dz."""
    if a.model == "linear":
        return d0 - dz
    # Gaussian-profile depression under linear diffusion, sigma0 = D/4:
    # depth = d0 / (1 + 2 kappa t / sigma0^2) = d0 / (1 + 32 kappa t / D^2)
    return d0 / (1.0 + 32.0 * dz / (D * 1000.0) ** 2)


def state(s, ns, ee, el, unit="N", tr=None, h=0.0):
    """Alive mask and current depth for snapshot ns (craters 0..ns-1 formed)."""
    n, area = s["n"], s["area"]
    idx = np.arange(n)
    T = np.full(n, -1.0)
    T[:ns] = t_of_n1((ns - idx[:ns]) / area) * 1000.0  # formation time, Myr before present
    Z = np.zeros(n)
    Z[:ns] = zcum(T[:ns], ee, el)
    d0 = s["d0"]
    alive = (idx < ns) & (s["imp"] >= ns)
    # ejecta burial, using the older crater's depth when the later crater formed
    i, j, th = s["ej_i"], s["ej_j"], s["ej_t"]
    k = j < ns
    i, j, th = i[k], j[k], th[k]
    depth_then = depth_after(d0[i], s["D"][i], Z[i] - Z[j])
    bur = th >= depth_then
    alive[np.unique(i[bur])] = False
    fill = np.zeros(n)
    if unit == "H":
        pre = (T > tr * 1000.0) & (idx < ns)
        rem = depth_after(d0, s["D"], Z - zcum(np.full(n, tr * 1000.0), ee, el))
        alive[pre & (rem <= h)] = False
        fill[pre] = h
    cur = depth_after(d0, s["D"], Z) - fill
    alive &= cur > a.qmin * d0
    if a.thin12:  # age-blind loss of 1-2 km craters (catalog incompleteness), matched to observed
        small = alive & (s["D"] < 2)
        target = obs[unit]["sfd"][0] * s["area"] / AREA[unit]
        if small.sum() > target and a.thin12 == 2:  # oldest first instead of at random
            idx_s = np.nonzero(small)[0]
            alive[idx_s[: len(idx_s) - int(target)]] = False
        elif small.sum() > target:
            alive &= ~(small & (s["u"] > target / small.sum()))
    return alive, cur


def match_ns(s, ee, el, target, unit="N", tr=None, h=0.0):
    """Snapshot whose surviving craters >= 16 km match the observed count."""
    big = s["D"] >= 16
    scale = AREA[unit] / s["area"]
    lo, hi = int(0.2 * s["n"]), s["n"]
    for _ in range(11):
        mid = (lo + hi) // 2
        al, _ = state(s, mid, ee, el, unit, tr, h)
        if (al & big).sum() * scale < target:
            lo = mid
        else:
            hi = mid
    return hi


def evaluate(ee, el, unit="N", tr=None, h=0.0, ns_fixed=None, use=None):
    sfd = np.zeros(len(BINS) - 1); qh = np.zeros(len(QB) - 1); out = []
    use = sims if use is None else use
    for s in use:
        ns = ns_fixed[s["name"]] if ns_fixed else match_ns(s, ee, el, obs["N"]["N16"])
        al, cur = state(s, ns, ee, el, unit, tr, h)
        Da = s["D"][al]
        sfd += np.histogram(Da, BINS)[0] * AREA[unit] / s["area"]
        sel = al & (s["D"] >= 5)
        qh += np.histogram(np.clip(cur[sel] / s["d0env"][sel], 0, 9.99), QB)[0]
        out.append(ns)
    sfd /= len(use)
    o = obs[unit]
    fb = BINS[:-1] >= a.fit_min_km
    return {"dev_sfd": poisson_dev(o["sfd"][fb], sfd[fb]), "dev_q": multinom_dev(o["qhist"], qh),
            "sfd": sfd.tolist(), "qfrac": (qh / qh.sum()).tolist(),
            "ns": {s["name"]: v for s, v in zip(use, out)}}


def depth_from(o, y, n):
    d = np.zeros(n, dtype=np.int64)
    for _ in range(64):
        new = d.copy()
        np.maximum.at(new, y, d[o] + 1)
        if np.array_equal(new, d):
            break
        d = new
    return d


def deg_decisions(s, cur):
    """Overlap decisions drawn from the observed class-dependent rates."""
    rng = np.random.default_rng(21)
    q = np.clip(cur / s["d0env"], 0, 9.99)
    qk = np.clip(np.searchsorted(QC, q, side="right") - 1, 0, len(QC) - 2)
    cum = np.cumsum(p_deg_q[qk], axis=1)
    deg = (rng.random(len(q))[:, None] > cum).sum(1)  # 0..3 = class 1..4
    o, y = s["older"], s["younger"]
    big = np.where(s["D"][o] >= s["D"][y], o, y); sml = np.where(big == o, y, o)
    p = dec_tab[deg[big], deg[sml]]
    small = s["D"][sml] < 3
    p = np.where(small, np.interp(np.log(s["D"][sml]), np.log([1.4, 2.8]), [0.69, 0.59]), p)
    return s["contain"] | (rng.random(len(o)) < p)


def graph_stats(ee, el, ns_map, unit="N", tr=None, h=0.0):
    res = []
    for s in sims:
        al, cur = state(s, ns_map[s["name"]], ee, el, unit, tr, h)
        row = {}
        dmask = deg_decisions(s, cur) if a.decision == "deg" else s["dec"]
        for name, emask in (("containment", s["contain"]), ("contain+decided_overlap", dmask)):
            m = al[s["older"]] & al[s["younger"]] & emask
            o, y = s["older"][m], s["younger"][m]
            d = depth_from(o, y, s["n"])
            indeg = np.bincount(y, minlength=s["n"]); outdeg = np.bincount(o, minlength=s["n"])
            linked = al & ((indeg + outdeg) > 0)
            sc = AREA[unit] / s["area"]
            row[name] = {"depth_counts": (np.bincount(d[linked], minlength=8)[:8] * sc).round(1).tolist(),
                         "indeg_counts": (np.bincount(indeg[linked], minlength=8)[:8] * sc).round(1).tolist()}
        res.append(row)
    return res


t0 = time.time()
if a.model == "linear":
    EE = np.geomspace(0.5, 30, 10)     # m/Myr before T_b
    EL = np.geomspace(1e-4, 0.1, 10)   # m/Myr after T_b
else:
    EE = np.geomspace(1e2, 3e5, 10)    # m^2/Myr before T_b
    EL = np.geomspace(0.1, 300, 10)    # m^2/Myr after T_b
if a.ee:
    lo_, hi_, k_ = a.ee.split(","); EE = np.geomspace(float(lo_), float(hi_), int(k_))
if a.el:
    lo_, hi_, k_ = a.el.split(","); EL = np.geomspace(float(lo_), float(hi_), int(k_))
grid = []
if a.best:
    EE, EL = [float(a.best.split(",")[0])], [float(a.best.split(",")[1])]
for ee in EE:
    for el in EL:
        r = evaluate(ee, el, use=sims[:a.grid_sims])
        grid.append({"eps_early": ee, "eps_late": el, **{k: r[k] for k in ("dev_sfd", "dev_q", "ns")}})
    print(f"eps_early {ee:.2f} done, {time.time() - t0:.0f} s")
G = pd.DataFrame(grid)
G["dev"] = G.dev_sfd + G.dev_q
best = G.loc[G.dev.idxmin()]
print(G.sort_values("dev").head(8)[["eps_early", "eps_late", "dev_sfd", "dev_q", "dev"]].round(2))
ee, el = float(best.eps_early), float(best.eps_late)
fitN = evaluate(ee, el)
nsN = fitN["ns"]
print("Noachian best", ee, el, "tau0", {k: v / sims[0]["area"] for k, v in nsN.items()})

# Hesperian units: same substrate and rates, resurfaced at T_r with fill h
hgrid = []
for tr in (3.5, 3.6, 3.7):
    for h in (0.05, 0.1, 0.2, 0.4, 0.8, 1.6):
        r = evaluate(ee, el, "H", tr, h, ns_fixed=nsN)
        hgrid.append({"t_r": tr, "h_km": h, "dev_sfd": r["dev_sfd"], "dev_q": r["dev_q"]})
H = pd.DataFrame(hgrid); H["dev"] = H.dev_sfd + H.dev_q
hb = H.loc[H.dev.idxmin()]
print(H.sort_values("dev").head(5).round(2))
fitH = evaluate(ee, el, "H", float(hb.t_r), float(hb.h_km), ns_fixed=nsN)

out = {"tb_Ga": a.tb, "qmin": a.qmin, "grid": G.to_dict("records"), "best": {"eps_early": ee, "eps_late": el},
       "fit_N": fitN, "obs_N": {k: v.tolist() if hasattr(v, "tolist") else v for k, v in obs["N"].items()},
       "graph_N": graph_stats(ee, el, nsN),
       "hgrid": H.to_dict("records"), "best_H": {"t_r": float(hb.t_r), "h_km": float(hb.h_km)},
       "fit_H": fitH, "obs_H": {k: v.tolist() if hasattr(v, "tolist") else v for k, v in obs["H"].items()},
       "graph_H": graph_stats(ee, el, nsN, "H", float(hb.t_r), float(hb.h_km)),
       "bins_km": BINS.tolist(), "q_bins": QB.tolist()}
out["model"] = a.model; out["thin12"] = a.thin12; out["decision"] = a.decision
json.dump(out, open(os.path.join(a.out, f"sim_erosion_fit_{a.tag}.json"), "w"), indent=1, default=float)
print("done", round(time.time() - t0), "s")
