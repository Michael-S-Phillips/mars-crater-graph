"""Method figure: the rim-continuity test on one real overlapping pair."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mcgraph.dem import MolaDEM, PPD
from mcgraph.evidence import rim_profiles
from mcgraph.catalog import R_MARS_KM

# python fig_method.py WORK OUT.png [BIG_ID SMALL_ID] [THEMIS_CUTOUT.npz]
work, out_png = sys.argv[1], sys.argv[2]
themis_path = next((x for x in sys.argv[3:] if x.endswith(".npz")), None)
ids = [x for x in sys.argv[3:] if not x.endswith(".npz")]
big_id, small_id = (ids[0], ids[1]) if len(ids) > 1 else ("07-0-000580", "07-0-000576")
df = pd.read_parquet(os.path.join(work, "catalog.parquet"))
dem = MolaDEM(os.path.join(work, "mola128_global.i2"))
b = int(df.index[df.CRATER_ID == big_id][0]); s = int(df.index[df.CRATER_ID == small_id][0])
prof = rim_profiles(df, [b, s], dem)
OLD, YOUNG, INK2, MUTED = "#e0612e", "#2a78d6", "#52514e", "#898781"
kpd = R_MARS_KM * np.pi / 180

def inside(k, m, thr):
    la1, lo1 = np.radians(prof[k]["lat"]), np.radians(prof[k]["lon"])
    la2, lo2 = np.radians(df.lat.iat[m]), np.radians(df.lon.iat[m])
    d = np.arccos(np.clip(np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2), -1, 1)) * R_MARS_KM
    return d < thr * df.r_km.iat[m], d > 1.25 * df.r_km.iat[m]

clat = 0.5 * (df.lat.iat[b] + df.lat.iat[s]); clon = 0.5 * (df.lon.iat[b] + df.lon.iat[s])
cosl = np.cos(np.radians(clat))
hw = 1.15 * (df.r_km.iat[b] + df.r_km.iat[s] + 0.5 * 0)  # km
n = 600 if themis_path else 360
ys = np.linspace(hw, -hw, n); xs = np.linspace(-hw, hw, n)
X, Y = np.meshgrid(xs, ys)
z = dem.sample(clat + Y / kpd, clon + X / (kpd * cosl))
gy, gx = np.gradient(z, (2 * hw / n) * 1000)
gy = -gy
slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
A, E = np.radians(315), np.radians(40)
hs = np.sin(E) * np.cos(slope) + np.cos(E) * np.sin(slope) * np.cos(A - aspect)
if themis_path:
    from mcgraph.basemap import load_cutout, sample_cutout, blend
    th = sample_cutout(load_cutout(themis_path), clat + Y / kpd, clon + X / (kpd * cosl))
    base, vmin, vmax = blend(th, hs), 0.0, 1.0
    base_name = "THEMIS day IR + MOLA relief"
else:
    base, vmin, vmax = hs, -0.2, 1.1
    base_name = "MOLA shaded relief"

fig = plt.figure(figsize=(10, 4.4), dpi=170)
ax = fig.add_axes([0.04, 0.1, 0.38, 0.82])
ax.imshow(base, cmap="gray", extent=(-hw, hw, -hw, hw), vmin=vmin, vmax=vmax)
for k, m, col, lab in ((b, s, OLD, "older (breached)"), (s, b, YOUNG, "younger (intact)")):
    xk = (((prof[k]["lon"] - clon + 180) % 360) - 180) * kpd * cosl; yk = (prof[k]["lat"] - clat) * kpd
    ax.plot(np.r_[xk, xk[0]], np.r_[yk, yk[0]], color=col, lw=0.9, ls=(0, (3, 2)))
    ins, _ = inside(k, m, 0.85)
    xi = np.where(ins, xk, np.nan); yi = np.where(ins, yk, np.nan)
    ax.plot(xi, yi, color=col, lw=3.2, solid_capstyle="butt")
    cx = (((df.lon.iat[k] - clon + 180) % 360) - 180) * kpd * cosl; cy = (df.lat.iat[k] - clat) * kpd
    ax.text(cx, cy, f"{df.CRATER_ID.iat[k]}\n{df.D_km.iat[k]:.0f} km, {lab}", color="white", fontsize=7,
            ha="center", va="center", bbox=dict(boxstyle="round,pad=0.25", fc=col, ec="none", alpha=0.9))
ax.set_xlim(-hw, hw); ax.set_ylim(-hw, hw); ax.set_aspect("equal")
ax.set_xlabel("km", fontsize=8, color=INK2); ax.tick_params(labelsize=7, colors=MUTED)
ax.set_title(f"(a) {base_name}; thick arcs = rim inside the partner", fontsize=8.5, loc="left", color="#0b0b0b")

ax2 = fig.add_axes([0.52, 0.13, 0.45, 0.76])
for k, m, col, lab in ((b, s, OLD, "older"), (s, b, YOUNG, "younger")):
    az = np.degrees(prof[k]["az"]); pr = prof[k]["prom"]
    ins, ref = inside(k, m, 0.85)
    P = np.median(pr[ins]) / np.median(pr[ref])
    ax2.plot(az, pr, color=col, lw=1.0, alpha=0.9, label=f"{df.CRATER_ID.iat[k]} ({lab}): P = {P:.2f}")
    ax2.scatter(az[ins], pr[ins], s=7, color=col, zorder=3)
ax2.axhline(0, color=MUTED, lw=0.8)
ax2.set_xlim(0, 360); ax2.set_xticks(range(0, 361, 90))
ax2.set_xlabel("Azimuth from crater center (degrees from north)", fontsize=8, color=INK2)
ax2.set_ylabel("Rim prominence (m)", fontsize=8, color=INK2)
ax2.tick_params(labelsize=7, colors=MUTED)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)
ax2.legend(fontsize=7, frameon=False, loc="lower right")
ax2.set_title("(b) Rim prominence; dots = azimuths on the intrusion arc", fontsize=8.5, loc="left", color="#0b0b0b")
fig.savefig(out_png)
print("ok", out_png)
