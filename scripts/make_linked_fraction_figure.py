"""Figure 3 as a static image: fraction of cataloged craters linked into the default graph,
by diameter and by host-unit period (values as in the paper's Figure 3 table).

  python scripts/make_linked_fraction_figure.py --nodes outputs/nodes.csv.gz --fig figures
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--nodes", required=True); ap.add_argument("--fig", required=True)
a = ap.parse_args()

nodes = pd.read_csv(a.nodes, low_memory=False, usecols=["D_km", "unit_period", "in_graph"])
edges = [1, 2, 5, 10, 20, 50, np.inf]
labels = ["1–2", "2–5", "5–10", "10–20", "20–50", "≥50"]
nodes["bin"] = pd.cut(nodes.D_km, edges, right=False, labels=labels)
frac = (nodes.groupby(["unit_period", "bin"], observed=False).in_graph.mean() * 100).unstack()
print(frac.round(1))

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SERIES = [("N", "Noachian units", "#2a78d6", "o"), ("H", "Hesperian units", "#eb6834", "s"),
          ("A", "Amazonian units", "#1baf7a", "^")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.8, 4.0))
x = np.arange(len(labels))
for key, name, col, mk in SERIES:
    y = frac.loc[key, labels].to_numpy(float)
    ax.plot(x, y, color=col, lw=2, marker=mk, ms=7, mec="white", mew=1.2, zorder=3)
    dy = {"N": 3.5, "H": -3.5, "A": 0}[key]  # keep the two top labels apart
    ax.text(x[-1] + 0.2, y[-1] + dy, f"{name} ({y[-1]:.1f}% at ≥50 km)", color=INK2, fontsize=8.5, va="center")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_xlim(-0.3, len(labels) + 1.9)
ax.set_ylim(0, 108); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlabel("Crater diameter (km)"); ax.set_ylabel("Craters linked into the graph (%)")
ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
fig.tight_layout()
out = os.path.join(a.fig, "fig_linked_fraction.png")
fig.savefig(out, dpi=200)
print("written", out)
