"""Blinded CTX validation sheets: raw CTX + CTX with circles labeled A/B (random)."""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chips_dir, sample_pq, catalog_pq, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)
s = pd.read_parquet(sample_pq).sample(frac=1, random_state=5)
df = pd.read_parquet(catalog_pq)
rng = np.random.default_rng(99)
K = 3389.5 * np.pi / 180
key = []; per = 8
s = s[[os.path.exists(os.path.join(chips_dir, f"{int(p)}.npz")) for p in s.pair_id]]
for sheet, start in enumerate(range(0, len(s), per)):
    rows = s.iloc[start:start + per]
    fig, axes = plt.subplots(4, 4, figsize=(12.8, 12.8), dpi=100)
    for ax in axes.ravel():
        ax.axis("off")
    for k, (_, r) in enumerate(rows.iterrows()):
        z = np.load(os.path.join(chips_dir, f"{int(r.pair_id)}.npz"))
        clat, clon, hw = float(z["clat"]), float(z["clon"]), float(z["hw_km"])
        img = z["ctx"].astype(float)
        v = img[img > 0]
        lo, hi = (np.percentile(v, [1, 99]) if v.size else (0, 255))
        labs = ["A", "B"]; rng.shuffle(labs)
        lab = {"big": labs[0], "small": labs[1]}
        key.append(dict(pair_id=int(r.pair_id), sheet=sheet, slot=k, big=lab["big"], small=lab["small"]))
        row, col = divmod(k, 2)
        for j in range(2):
            ax = axes[row, col * 2 + j]
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, extent=(-hw, hw, -hw, hw), origin="upper")
            if j == 1:
                for tag in ("big", "small"):
                    node = int(r[tag])
                    dlon = ((df.lon.iat[node] - clon + 180) % 360) - 180
                    x = dlon * np.cos(np.radians(clat)) * K; y = (df.lat.iat[node] - clat) * K
                    rr = df.r_km.iat[node]; t = np.linspace(0, 2 * np.pi, 500)
                    c = "#00e5ff" if lab[tag] == "A" else "#ff4fd8"
                    ax.plot(x + rr * np.cos(t), y + rr * np.sin(t), ls=(0, (4, 4)), lw=1.0, color=c, alpha=0.8)
                    ang = np.arctan2(-y, -x) if (abs(x) + abs(y)) > 1e-6 else np.pi / 4
                    lx = np.clip(x + rr * np.cos(ang) * 0.6, -hw * 0.88, hw * 0.88)
                    ly = np.clip(y + rr * np.sin(ang) * 0.6, -hw * 0.88, hw * 0.88)
                    ax.text(lx, ly, lab[tag], color=c, fontsize=13, weight="bold", ha="center", va="center",
                            bbox=dict(boxstyle="round,pad=0.1", fc="black", alpha=0.55, lw=0))
            ax.set_xlim(-hw, hw); ax.set_ylim(-hw, hw)
            ax.set_title(f"#{sheet}.{k}  CTX  {2*hw:.1f} km" + ("" if j == 0 else "  circles"), fontsize=9)
    plt.tight_layout(); fig.savefig(os.path.join(out_dir, f"ctx_sheet_{sheet:02d}.png")); plt.close(fig)
pd.DataFrame(key).to_csv(os.path.join(out_dir, "KEY_do_not_open.csv"), index=False)
print("sheets", sheet + 1, "pairs", len(key))
