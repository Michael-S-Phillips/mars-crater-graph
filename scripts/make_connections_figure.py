"""Figure: number of craters with at least N relations (Figure 5 of the paper).

  python scripts/make_connections_figure.py --out outputs --fig figures
"""
import argparse, os
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--fig", required=True)
a = ap.parse_args()
root = a.out.rstrip("/") + "/"
e = pd.read_csv(root + "edges.csv"); e = e[e.kept]
n = pd.read_csv(root + "nodes.csv.gz", low_memory=False)
deg = pd.concat([e.older, e.younger]).value_counts()
vals = np.sort(deg.to_numpy())[::-1]                     # most-connected first
N = np.unique(vals)[::-1]                                # distinct N, high to low
cum = np.array([(vals >= k).sum() for k in N])           # craters with >= N connections
INK, INK2, MUTED, GRID, BLUE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=180)
ax.step(N, cum, where="post", color=BLUE, lw=2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
ax.set_xlim(N[0] * 1.6, 0.8); ax.set_ylim(0.7, 2.5e5)
ax.grid(color=GRID, lw=0.6, which="major"); ax.set_axisbelow(True)
ax.set_xlabel("Number of connections, N (relations to other craters; most to fewest)")
ax.set_ylabel("Craters with at least N connections")
ax.set_title(f"{len(deg):,} linked craters, {len(e):,} relations", loc="left", fontsize=10, color=INK)
# label the most-connected craters
for rank, (node, k) in enumerate(deg.head(4).items()):
    r = n.loc[node]
    name = "Ladon basin" if r.CRATER_ID == "07-4-043926" else r.CRATER_ID
    ax.annotate(f"{name} ({r.D_km:.0f} km): {k:,}", (k, rank + 1), xytext=(5200, 40 * 3.2 ** (3 - rank)),
                textcoords="data", fontsize=7.5, color=INK2, va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=2, shrinkB=3))
    ax.plot(k, rank + 1, "o", ms=3.5, color=BLUE)
for k in (10, 1):
    c = (vals >= k).sum()
    ax.annotate(f"N ≥ {k}: {c:,} craters", (k, c), xytext=(-8 if k == 1 else 8, -14 if k == 1 else -16),
                textcoords="offset points", fontsize=7.5, color=INK2, ha="right" if k == 1 else "left")
    ax.plot(k, c, "o", ms=3.5, color=BLUE)
ax.text(0.99, 0.02, f"{(n.in_graph == False).sum():,} cataloged craters have N = 0 (not shown)",
        transform=ax.transAxes, ha="right", fontsize=7, color=MUTED)
fig.tight_layout()
fig.savefig(os.path.join(a.fig, "fig_connections_cumulative.png"))
print(len(N), N[:5], cum[:5], cum[-1])
