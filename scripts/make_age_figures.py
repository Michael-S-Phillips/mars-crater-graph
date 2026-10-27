"""Figures for the age-dating analysis.

  python scripts/make_age_figures.py --out outputs --fig figures
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.ages import t_of_n1, epoch_boundaries_tau
from mcgraph.chronology import log_lik

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--fig", required=True)
a = ap.parse_args()

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE, ORANGE = "#2a78d6", "#eb6834"
DEG_COLORS = {1: "#0d366b", 2: "#1c5cab", 3: "#3987e5", 4: "#86b6ef"}  # darker = more degraded
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})

r = pd.read_parquet(os.path.join(a.out, "interior_counts.parquet"))
n = pd.read_csv(os.path.join(a.out, "nodes.csv.gz"), low_memory=False)
r = r.join(n[["DEG_RIM"]])
test = json.load(open(os.path.join(a.out, "age_pair_test.json")))
eb = epoch_boundaries_tau()

# ---------------- Figure 4: interior densities, degradation, and order ----------------
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw={"width_ratios": [1.25, 1, 1]})
ax = axs[0]
dms = ["1", "1.25", "1.5", "2", "3", "4", "5"]
x = np.array([float(d) for d in dms])
rng = np.random.default_rng(0)
labels = {1: "1 (most degraded)", 2: "2", 3: "3", 4: "4 (freshest)"}
deg_rows = []
for deg in [1, 2, 3, 4]:
    g = r[r.DEG_RIM == deg]
    K = np.stack([g[f"k_{d}"].to_numpy() for d in dms], 1)
    L = np.stack([g[f"L_{d}"].to_numpy() for d in dms], 1)
    t = K.sum(0) / L.sum(0)
    bs = []
    for _ in range(300):
        i = rng.integers(0, len(g), len(g))
        bs.append(K[i].sum(0) / L[i].sum(0))
    lo, hi = np.quantile(bs, [0.025, 0.975], axis=0)
    deg_rows += [{"deg": deg, "dmin_km": float(d), "n": len(g), "k": int(K[:, j].sum()),
                  "tau": t[j], "tau_lo95": lo[j], "tau_hi95": hi[j], "T_Ga": float(t_of_n1(t[j])),
                  "T_lo95": float(t_of_n1(lo[j])), "T_hi95": float(t_of_n1(hi[j]))} for j, d in enumerate(dms)]
    c = DEG_COLORS[deg]
    ax.fill_between(x, lo * 1e3, hi * 1e3, color=c, alpha=0.18, lw=0)
    ax.plot(x, t * 1e3, color=c, lw=2, marker="o", ms=4, label=f"{labels[deg]}  (n={len(g):,})")
pd.DataFrame(deg_rows).to_csv(os.path.join(a.out, "degradation_stacked_ages.csv"), index=False)
for name in ["EA/LH", "EH/LN", "LN/MN"]:
    yv = eb[name] * 1e3
    ax.axhline(yv, color=GRID, lw=0.8, zorder=0)
    ax.text(5.15, yv, f"{name}\n{t_of_n1(eb[name]):.2f} Ga", fontsize=6.5, color=MUTED, va="center")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([1, 1.5, 2, 3, 4, 5]); ax.set_xticklabels(["1", "1.5", "2", "3", "4", "5"])
ax.minorticks_off()
ax.set_yticks([0.2, 0.5, 1, 2, 5, 10, 20]); ax.set_yticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
ax.set_ylim(0.2, 20); ax.set_xlim(0.95, 6.6)
ax.set_xlabel("Smallest crater counted, D$_{min}$ (km)")
ax.set_ylabel("Interior crater density, N(1)-equivalent\nper 1,000 km$^2$")
ax.legend(title="Rim degradation class", fontsize=7, title_fontsize=7, loc="lower right",
          frameon=False, bbox_to_anchor=(1.0, 0.0), ncol=2, columnspacing=0.8, handlelength=1.5)
ax.set_title("a  Interior densities by rim degradation", loc="left", fontsize=9.5)

ax = axs[1]
bd = test["by_dmin"]
xs = np.array([float(k) for k in bd])
rat = np.array([v["ratio"] for v in bd.values()])
rlo = np.array([v["ratio_95"][0] for v in bd.values()]); rhi = np.array([v["ratio_95"][1] for v in bd.values()])
ax.axhline(1, color=MUTED, lw=1, ls=(0, (3, 3)))
ax.text(3.05, 1.04, "no difference", fontsize=7, color=MUTED, ha="right", va="bottom")
ax.errorbar(xs, rat, yerr=[rat - rlo, rhi - rat], fmt="o", color=BLUE, ms=6, lw=1.5, capsize=0)
for xx, yy in zip(xs, rat):
    ax.text(xx * 1.06, yy, f"{yy:.2f}", fontsize=7.5, color=INK2, va="center")
ax.set_xscale("log"); ax.set_xticks(xs); ax.set_xticklabels([f"{v:g}" for v in xs]); ax.minorticks_off()
ax.set_xlim(0.85, 3.8); ax.set_ylim(0.8, 2.7)
ax.set_xlabel("D$_{min}$ (km)")
ax.set_ylabel("Interior density, older ÷ younger crater")
ax.set_title("b  Older crater's interior is denser", loc="left", fontsize=9.5)

ax = axs[2]
fi = np.array([v["frac_inverted"] for v in bd.values()]) * 100
flo = np.array([v["frac_inverted_95"][0] for v in bd.values()]) * 100
fhi = np.array([v["frac_inverted_95"][1] for v in bd.values()]) * 100
ax.errorbar(xs, fi, yerr=[fi - flo, fhi - fi], fmt="o", color=ORANGE, ms=6, lw=1.5, capsize=0)
for xx, yy in zip(xs, fi):
    ax.text(xx * 1.06, yy, f"{yy:.1f}%", fontsize=7.5, color=INK2, va="center")
ax.set_xscale("log"); ax.set_xticks(xs); ax.set_xticklabels([f"{v:g}" for v in xs]); ax.minorticks_off()
ax.set_xlim(0.85, 3.8); ax.set_ylim(-0.5, 17)
ax.set_xlabel("D$_{min}$ (km)")
ax.set_ylabel("Edges whose older crater has the\nsparser interior (%)")
ax.set_title("c  Order inverted by the counts", loc="left", fontsize=9.5)
for axx in axs:
    axx.grid(axis="y", color=GRID, lw=0.6); axx.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(a.fig, "fig4_interior_ages.png"), dpi=200)
plt.close(fig)

# ---------------- Figure 5: a superposition sequence (Copernicus) ----------------
o = pd.read_parquet(os.path.join(a.out, "crater_ages_main.parquet"))
z = np.load(os.path.join(a.out, "crater_ages_main_samples.npz"))
grid, snodes, S = z["grid"], z["nodes"], z["idx"]
pr = json.load(open(os.path.join(a.out, "crater_ages_main_prior.json")))
g = np.array(next(iter(pr["priors"].values())))
host = o.index[o.CRATER_ID == "04-1-001446"][0]
e = pd.read_csv(os.path.join(a.out, "edges.csv")); e = e[e.kept]
kids = [y for y in e.younger[e.older == host] if y in o.index]
rows = [host] + sorted(kids, key=lambda v: -o.D_km[v])
spos = {v: i for i, v in enumerate(snodes)}
XMIN = 2.0
fig, ax = plt.subplots(figsize=(7.6, 3.6))
for j, v in enumerate(rows):
    yv = len(rows) - 1 - j
    ll = log_lik([o.k[v]], [o.Lambda[v]], grid)[0] + np.log(g + 1e-300)
    p = np.exp(ll - ll.max()); p /= p.sum(); c = np.cumsum(p)
    qu = t_of_n1(grid[[np.searchsorted(c, q) for q in (0.025, 0.16, 0.5, 0.84, 0.975)]])
    ts = t_of_n1(grid[S[:, spos[v]]])
    qc = np.quantile(ts, [0.025, 0.16, 0.5, 0.84, 0.975])
    for q, off, col in ((qu, 0.17, MUTED), (qc, -0.17, BLUE)):
        if q[0] < XMIN:
            ax.plot(XMIN + 0.02, yv + off, marker="<", color=col, ms=4, mew=0, clip_on=False)
        ax.plot([max(q[0], XMIN + 0.02), q[4]], [yv + off] * 2, color=col, lw=1, solid_capstyle="round")
        ax.plot([max(q[1], XMIN + 0.02), q[3]], [yv + off] * 2, color=col, lw=4, solid_capstyle="round")
        ax.plot(q[2], yv + off, "o", color="white", ms=3.5, mec=col, mew=1.2)
    nm = "Copernicus" if v == host else o.CRATER_ID[v]
    ax.text(XMIN - 0.03, yv, f"{nm}  {o.D_km[v]:.0f} km, k = {o.k[v]}", fontsize=7.5,
            color=INK if v == host else INK2, va="center", ha="right")
ax.axhline(len(rows) - 1.5, color=GRID, lw=0.8)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlim(XMIN, 3.95)
ax.set_xlabel("Model formation age (Ga), Hartmann & Neukum (2001) chronology")
ax.plot([], [], color=MUTED, lw=4, label="interior count only")
ax.plot([], [], color=BLUE, lw=4, label="count + superposition order")
ax.legend(loc="upper left", fontsize=7.5, frameon=False, bbox_to_anchor=(0.0, 1.02))

fig.suptitle("Copernicus crater (49°S, 169°W) and the seven craters ≥10 km superposed on it",
             x=0.02, ha="left", fontsize=9.5)
fig.tight_layout()
fig.subplots_adjust(left=0.30, top=0.9)
fig.savefig(os.path.join(a.fig, "fig5_sequence_example.png"), dpi=200)
print("written")
