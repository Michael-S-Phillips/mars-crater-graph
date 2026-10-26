"""Figures: global overview, a regional Hasse-diagram zoom, and validation summary.

  python scripts/make_figures.py --work WORK --out outputs --fig figures
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.dem import MolaDEM, PPD

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--fig", required=True)
ap.add_argument("--region", default="-37.5,-33.0,358.5,3.5",
                help="lat_min,lat_max,lon_min,lon_max for the zoom panel (lon may wrap)")
a = ap.parse_args()
os.makedirs(a.fig, exist_ok=True)

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
DEPTH_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # ordinal blue ramp
OLDER, YOUNGER = "#eb6834", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED})

nodes = pd.read_csv(os.path.join(a.out, "nodes.csv.gz"), index_col="node", low_memory=False)
edges = pd.read_csv(os.path.join(a.out, "edges.csv"))
edges = edges[edges.kept & edges.in_hasse]
dem = MolaDEM(os.path.join(a.work, "mola128_global.i2"))


def hillshade(z, dx_m, dy_m, az=315, alt=40):
    gy, gx = np.gradient(z.astype(float), dy_m, dx_m)
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
    A, E = np.radians(az), np.radians(alt)
    return np.sin(E) * np.cos(slope) + np.cos(E) * np.sin(slope) * np.cos(A - aspect)


# ---------------- Figure 1: global overview ----------------
step = 16  # 8 ppd
z = np.asarray(dem.a[::step, ::step], dtype=float)
px_m = 3396000 * np.pi / 180 / (PPD / step)
hs = hillshade(z, px_m, px_m)
g = nodes[nodes.in_graph].copy()
fig, ax = plt.subplots(figsize=(12, 6.4), dpi=150)
ax.imshow(hs, cmap="gray", extent=(0, 360, -88, 88), vmin=-0.2, vmax=1.2, alpha=0.55,
          interpolation="bilinear")
lon = g.lon.to_numpy(); lat = g.lat.to_numpy()
dep = g.depth.astype(int).clip(upper=3).to_numpy()
size = np.clip(g.D_km.to_numpy() / 4, 0.4, 60)
order = np.argsort(dep)
ax.scatter(lon[order], lat[order], s=size[order], c=np.array(DEPTH_COLORS)[dep[order]],
           linewidths=0, alpha=0.85)
ax.set_xlim(0, 360); ax.set_ylim(-88, 88)
ax.set_xticks(range(0, 361, 60)); ax.set_yticks(range(-60, 61, 30))
ax.set_xlabel("East longitude (°)"); ax.set_ylabel("Latitude (°)")
for s in ax.spines.values():
    s.set_visible(False)
counts = g.depth.astype(int).clip(upper=3).value_counts().sort_index()
labels = ["depth 0  (nothing older known beneath)", "depth 1", "depth 2", "depth 3+"]
handles = [Line2D([], [], marker="o", ls="", color=DEPTH_COLORS[i], markersize=6,
                  label=f"{labels[i]}  n={counts.get(i, 0):,}") for i in range(4)]
ax.legend(handles=handles, loc="lower left", frameon=True, framealpha=0.9, fontsize=8,
          title="Longest chain of older craters beneath", title_fontsize=8)
ax.set_title(f"{len(g):,} Mars craters linked by superposition "
             f"({int(edges.shape[0]):,} non-redundant older-to-younger relations)",
             loc="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(os.path.join(a.fig, "fig1_global_depth_map.png")); plt.close(fig)

# ---------------- Figure 2: regional Hasse diagram ----------------
la0, la1, lo0, lo1 = map(float, a.region.split(","))
if lo1 < lo0:
    lo1 += 360
clat, clon = 0.5 * (la0 + la1), 0.5 * (lo0 + lo1)
kpd = 3389.5 * np.pi / 180
cosl = np.cos(np.radians(clat))
def to_xy(lat, lon):
    dl = ((np.asarray(lon) - clon + 180) % 360) - 180
    return dl * cosl * kpd, (np.asarray(lat) - clat) * kpd
ny = int((la1 - la0) * PPD); nx_ = int((lo1 - lo0) * PPD)
LA, LO = np.meshgrid(np.linspace(la1, la0, ny), np.linspace(lo0, lo1, nx_), indexing="ij")
zr = dem.sample(LA, np.mod(LO, 360))
hs2 = hillshade(zr, kpd * 1000 / PPD * cosl, kpd * 1000 / PPD)
x0, y0 = to_xy(la0, lo0); x1, y1 = to_xy(la1, lo1)
nl = np.mod(nodes.lon.to_numpy(), 360)
inlon = ((nl - lo0) % 360) <= (lo1 - lo0)
sel = nodes[nodes.in_graph.to_numpy() & nodes.lat.between(la0, la1).to_numpy() & inlon]
ids = set(sel.index)
ee = edges[edges.older.isin(ids) & edges.younger.isin(ids)]
W = 10; H = W * (y1 - y0) / (x1 - x0)
fig, ax = plt.subplots(figsize=(W, H + 0.6), dpi=150)
ax.imshow(hs2, cmap="gray", extent=(x0, x1, y0, y1), vmin=-0.1, vmax=1.1, alpha=0.8)
X, Y = to_xy(sel.lat, sel.lon)
for (n, r), xc, yc in zip(sel.iterrows(), X, Y):
    d = int(min(r.depth, 3))
    ax.add_patch(Circle((xc, yc), r.D_km / 2, fill=False, lw=1.0, ec=DEPTH_COLORS[d]))
pos = {n: (xc, yc) for n, xc, yc in zip(sel.index, X, Y)}
for r in ee.itertuples():
    col = INK if r.method == "containment" else OLDER
    ax.add_patch(FancyArrowPatch(pos[r.older], pos[r.younger], arrowstyle="-|>", mutation_scale=7,
                                 lw=0.7, color=col, alpha=0.9, shrinkA=0, shrinkB=1.5))
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
ax.set_xlabel(f"km east of {clon % 360:.1f}°E"); ax.set_ylabel(f"km north of {clat:.1f}°")
handles = [Line2D([], [], color=INK, lw=1, label="containment (geometric)"),
           Line2D([], [], color=OLDER, lw=1, label="rim overlap (rim-arc / topography test)")]
handles += [Line2D([], [], marker="o", ls="", mfc="none", mec=DEPTH_COLORS[i], markersize=7,
                   label=f"depth {i}{'+' if i == 3 else ''}") for i in range(4)]
ax.legend(handles=handles, loc="lower left", fontsize=7.5, framealpha=0.9)
ax.set_title(f"Arrows run older → younger (transitive reduction): {len(sel)} craters, {len(ee)} relations",
             loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(a.fig, "fig2_regional_hasse.png")); plt.close(fig)

# ---------------- Figure 3: validation summary ----------------
v = json.load(open(os.path.join(a.out, "validation.json")))
q = v["deg_ceiling_from_containment"]["q"]
label = {"arc": "rim-arc edges vs. degradation state",
         "topo+arc": "rim-arc + topography edges vs. degradation state",
         "topo": "topography-only edges vs. degradation (low; excluded)"}
rows = []
for m in sorted(v["deg_check_by_method"], key=lambda m: ["topo+arc", "arc", "topo"].index(m["method"])):
    rows.append((label[m["method"]], m["deg_agreement"], m["n_with_deg"]))
rows.append(("rim-arc test vs. topography test (independent data)", v["topo_vs_arc"]["agreement"], v["topo_vs_arc"]["n"]))
cv = v.get("ctx_visual")
if cv:
    rows.append(("blind visual check on CTX images", cv["agreement"], cv["oriented_and_judgeable"]))
rows.append(("younger crater is the one with layered ejecta", v["layered_ejecta_check"]["younger_is_the_layered_one"],
             v["layered_ejecta_check"]["n"]))
fig, ax = plt.subplots(figsize=(8.5, 3.8), dpi=150)
yv = np.arange(len(rows))[::-1]
for y_ in yv:
    ax.axhline(y_, color=GRID, lw=0.8, zorder=0)
ax.scatter([r[1] * 100 for r in rows], yv, s=46, color="#2a78d6", zorder=3)
for y_, r in zip(yv, rows):
    ax.text(r[1] * 100 + 1.2, y_ + 0.18, f"{r[1]*100:.1f}%  (n={r[2]:,})", va="bottom", fontsize=7.5, color=INK2)
ax.axvline(q * 100, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
ax.text(q * 100 - 0.6, -0.95, f"degradation-check ceiling, set by containment edges ({q*100:.1f}%)",
        fontsize=7, color=INK2, ha="right", va="center")
ax.set_yticks(yv); ax.set_yticklabels([r[0] for r in rows], fontsize=8, color=INK2)
ax.set_xlim(50, 109); ax.set_ylim(-1.4, len(rows) - 0.4); ax.set_xlabel("Agreement (%)")
for s_ in ("top", "right", "left"):
    ax.spines[s_].set_visible(False)
ax.tick_params(axis="y", length=0)
fig.suptitle("Agreement of the inferred older-to-younger order with independent evidence", x=0.01, ha="left",
             fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(a.fig, "fig3_validation.png")); plt.close(fig)
print("figures written")
