"""Figure: degradation-model fit on Noachian terrain and its prediction of chain depth.

  python scripts/make_erosion_figure.py --out outputs --fig figures
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--fig", required=True)
ap.add_argument("--final", default="final_diffusion")
ap.add_argument("--linear", default="linear_arc0.3")
ap.add_argument("--half", default="diffusion", help="fit run that used the 50% rim-loss rule")
a = ap.parse_args()
J = lambda t: json.load(open(os.path.join(a.out, f"sim_erosion_fit_{t}.json")))
fin, lin, half = J(a.final), J(a.linear), J(a.half)

INK, INK2, MUTED, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.9))

# (a) size-frequency distribution, incremental per factor-of-2 bin
ax = axs[0]
b = np.array(fin["bins_km"]); c = np.sqrt(b[:-1] * b[1:])
ax.plot(c, fin["obs_N"]["sfd"], "s", color=INK, ms=5, label="observed, Noachian units", zorder=3)
ax.plot(c, fin["fit_N"]["sfd"], "-o", color=BLUE, ms=3.5, lw=2, label="diffusive degradation (best fit)")
ax.plot(c, lin["fit_N"]["sfd"], "--", color=ORANGE, lw=1.5, label="linear depth loss (best fit)")
ax.axvspan(1, 2, color=GRID, alpha=0.5, lw=0)
ax.text(1.05, 12, "1–2 km matched by\nage-blind thinning", fontsize=6.5, color=MUTED)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Crater diameter (km)"); ax.set_ylabel("Craters per factor-of-2 diameter bin")
ax.set_xticks([1, 4, 16, 64, 256]); ax.set_xticklabels(["1", "4", "16", "64", "256"]); ax.minorticks_off()
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax.set_title("a  Size-frequency distribution", loc="left", fontsize=9.5, color=INK)

# (b) relative depth distribution, craters >= 5 km
ax = axs[1]
qb = np.array(fin["q_bins"]); qb[-1] = 1.2
qc = 0.5 * (qb[:-1] + qb[1:])
o = np.array(fin["obs_N"]["qhist"]); o = o / o.sum()
w = np.diff(qb)
ax.bar(qc, o / w, width=w * 0.92, color=GRID, edgecolor="none", label="observed (MOLA depths)")
ax.plot(qc, np.array(fin["fit_N"]["qfrac"]) / w, "-o", color=BLUE, ms=3.5, lw=2, label="diffusive degradation")
ax.plot(qc, np.array(lin["fit_N"]["qfrac"]) / w, "--", color=ORANGE, lw=1.5, label="linear depth loss")
ax.set_xlabel("Depth ÷ fresh depth (craters ≥5 km)"); ax.set_ylabel("Probability density")
ax.set_xlim(0, 1.2)
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax.set_title("b  Relative crater depth", loc="left", fontsize=9.5, color=INK)

# (c) chain depth: predictions, not fitted
ax = axs[2]
K = 6; x = np.arange(K)
obs_all = [13009, 44936, 4416, 152, 2, 0]
for run, col, lab, ls in ((half, MUTED, "craters lost at 50% rim loss", (0, (4, 2))),
                          (fin, BLUE, "craters lost at 30% rim loss", "-")):
    R = np.array([g["contain+decided_overlap"]["depth_counts"][:K] if "contain+decided_overlap" in g
                  else g["contain+sized_overlap"]["depth_counts"][:K] for g in run["graph_N"]])
    m = R.mean(0); ok = m > 0.3
    ax.plot(x[ok], m[ok], ls=ls, color=col, lw=2, marker="o", ms=3.5, label=lab)
    ax.vlines(x[ok], np.maximum(R.min(0)[ok], 0.5), R.max(0)[ok], color=col, lw=1)
ob = np.array(obs_all, float); ok = ob > 0
ax.plot(x[ok], ob[ok], "s", color=INK, ms=5, label="observed")
ax.set_yscale("log"); ax.set_ylim(0.5, 1e5); ax.set_xticks(x)
ax.set_xlabel("Chain depth"); ax.set_ylabel("Linked craters")
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax.set_title("c  Chain depth (predicted, not fitted)", loc="left", fontsize=9.5, color=INK)
for axx in axs:
    axx.grid(axis="y", color=GRID, lw=0.6); axx.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(a.fig, "fig7_erosion_fit.png"), dpi=200)
print("written")
