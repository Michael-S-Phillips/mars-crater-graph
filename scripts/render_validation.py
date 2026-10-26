"""Render blinded validation sheets: THEMIS day IR + DEM hillshade with the two
crater circles labeled A/B in random order. The key is written separately."""
import json, sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chips_dir, sample_pq, catalog_pq, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)
s = pd.read_parquet(sample_pq)
df = pd.read_parquet(catalog_pq)
rng = np.random.default_rng(2026)
s = s.sample(frac=1, random_state=11)  # shuffle so strata are mixed
K = 3389.5 * np.pi / 180


def hillshade(z, az=315, alt=35, px_km=1.0):
    gy, gx = np.gradient(z, px_km * 1000)
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy)
    a, e = np.radians(az), np.radians(alt)
    return np.sin(e) * np.cos(slope) + np.cos(e) * np.sin(slope) * np.cos(a - aspect)


key = []
per = 8
for sheet, start in enumerate(range(0, len(s), per)):
    rows = s.iloc[start:start + per]
    fig, axes = plt.subplots(4, 4, figsize=(12.8, 12.8), dpi=100)
    for ax in axes.ravel():
        ax.axis("off")
    for k, (_, r) in enumerate(rows.iterrows()):
        z = np.load(os.path.join(chips_dir, f"{int(r.pair_id)}.npz"))
        clat, clon, hw = float(z["clat"]), float(z["clon"]), float(z["hw_km"])
        labs = ["A", "B"]; rng.shuffle(labs)
        lab = {"big": labs[0], "small": labs[1]}
        key.append(dict(pair_id=int(r.pair_id), sheet=sheet, slot=k, big=lab["big"], small=lab["small"]))
        th = z["themis"].astype(float); dem = z["dem"].astype(float)
        dem[dem <= -32000] = np.nan
        hs = hillshade(np.nan_to_num(dem, nan=np.nanmean(dem)), px_km=2 * hw / dem.shape[0])
        row, col = divmod(k, 2)
        for j, (img, cmap) in enumerate(((th, "gray"), (hs, "gray"))):
            ax = axes[row, col * 2 + j]
            lo, hi = np.percentile(img, [1, 99])
            ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi, extent=(-hw, hw, -hw, hw), origin="upper")
            for tag, node, color in (("big", int(r.big), "#00e5ff"), ("small", int(r.small), "#ffd400")):
                la, lo_ = df.lat.iat[node], df.lon.iat[node]
                dlon = ((lo_ - clon + 180) % 360) - 180
                x = dlon * np.cos(np.radians(clat)) * K; y = (la - clat) * K
                rr = df.r_km.iat[node]
                th_ = np.linspace(0, 2 * np.pi, 400)
                c = "#00e5ff" if lab[tag] == "A" else "#ff4fd8"
                ax.plot(x + rr * np.cos(th_), y + rr * np.sin(th_), ls=(0, (4, 4)), lw=0.9, color=c, alpha=0.75)
                # label at the circle point nearest the window, offset inward
                ang = np.arctan2(-y, -x) if (x or y) else np.pi / 4
                lx, ly = x + rr * np.cos(ang) * 0.75, y + rr * np.sin(ang) * 0.75
                lx, ly = np.clip(lx, -hw * 0.9, hw * 0.9), np.clip(ly, -hw * 0.9, hw * 0.9)
                ax.text(lx, ly, lab[tag], color=c, fontsize=13, weight="bold", ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.1", fc="black", alpha=0.5, lw=0))
            ax.set_xlim(-hw, hw); ax.set_ylim(-hw, hw)
            if j == 0:
                ax.set_title(f"#{sheet}.{k}  ({2*hw:.0f} km wide)  THEMIS", fontsize=9)
            else:
                ax.set_title(f"#{sheet}.{k}  HRSC-MOLA hillshade", fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f"sheet_{sheet:02d}.png"))
    plt.close(fig)
pd.DataFrame(key).to_csv(os.path.join(out_dir, "KEY_do_not_open.csv"), index=False)
print("sheets", sheet + 1)
