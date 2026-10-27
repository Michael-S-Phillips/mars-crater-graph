"""Figure: chain depth on Noachian terrain, observed and simulated.

  python scripts/make_sim_figure.py --out outputs --fig figures [--erasure arc]
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--fig", required=True)
ap.add_argument("--erasure", default="arc")
a = ap.parse_args()
j = json.load(open(os.path.join(a.out, f"sim_compare_{a.erasure}.json")))
obs, sim = j["observed_noachian"], j["simulated"]

INK, INK2, MUTED, GRID, BLUE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
K = 8


def runs(filt, edges):
    return np.array([s["depth_counts"][:K] for s in sim if s["filter"] == filt and s["edges"] == edges], float)


panels = [("containment", "containment", "containment", "a  Containment edges only"),
          ("contain+decided_overlap", "contain+all_overlap", "contain+sized_overlap",
           "b  Containment and rim-overlap edges")]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.3), sharey=True)
for ax, (okey, e_all, e_dec, title) in zip(axs, panels):
    x = np.arange(K)
    for filt, ekey, col, lab, ls in (
            ("none", e_all, MUTED, "simulated: impacts only, every relation recorded", (0, (4, 2))),
            ("thin_oldest", e_dec, BLUE, "simulated: + erosion of oldest craters + pipeline decision rates", "-")):
        R = runs(filt, ekey)
        m = R.mean(0)
        ok = m > 0
        ax.plot(x[ok], m[ok], color=col, lw=2, ls=ls, marker="o", ms=4, label=lab)
        lo, hi = R.min(0), R.max(0)
        ax.vlines(x[ok], np.maximum(lo[ok], 0.5), hi[ok], color=col, lw=1)
    o = np.array(obs[okey]["depth_counts"][:K], float)
    ok = o > 0
    ax.plot(x[ok], o[ok], color=INK, lw=0, marker="s", ms=6, label="observed, Noachian units")
    for xx, yy in zip(x[ok], o[ok]):
        ax.text(xx + 0.15, yy * 1.15, f"{int(yy):,}", fontsize=7, color=INK2)
    ax.set_yscale("log"); ax.set_ylim(0.5, 4e5); ax.set_xlim(-0.4, K - 0.5)
    ax.set_xticks(x); ax.set_xlabel("Depth (longest chain of older craters beneath)")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=9.5, color=INK)
axs[0].set_ylabel("Linked craters, Noachian-unit area")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.text(0.01, 0.01, "Craters 1–362 km on 6.8×10⁷ km² (the area of Noachian units); simulation snapshot matched to the "
         "observed density of craters ≥16 km; bars span three simulation runs.", fontsize=6.5, color=MUTED)
fig.tight_layout(rect=(0, 0.04, 1, 0.93))
fig.savefig(os.path.join(a.fig, "fig6_chain_depth_simulation.png"), dpi=200)
print("written")
