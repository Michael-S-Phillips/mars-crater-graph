"""Step 4: validate edge orientations against independent evidence.

  python scripts/validate.py --work WORK --out outputs [--ctx validation/ctx_results.csv]

1. Degradation-state check. Robbins's rim degradation class (DEG_RIM, 4 = fresh)
   is an imperfect age proxy, so we first measure how often it agrees with
   CONTAINMENT edges, whose direction is geometrically certain. That agreement q
   (among pairs whose classes differ) is the ceiling any correct method can
   reach. For a method with observed agreement A, the implied accuracy is
   a = (A - (1 - q)) / (2q - 1), assuming degradation errors are independent of
   the method.
2. Cross-agreement between the topographic (E1) and rim-arc (E2) tests, which use
   independent data (MOLA vs. THEMIS-based rim tracing).
3. Layered-ejecta check: when exactly one crater of an oriented overlapping pair
   has mapped layered ejecta, that crater should usually be the younger one.
4. Blind visual check on CTX (if a results table is supplied).
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--ctx", default=None)
a = ap.parse_args()

df = pd.read_parquet(os.path.join(a.work, "catalog.parquet"))
p = pd.read_parquet(os.path.join(a.work, "pairs_evidence.parquet"))
deg = df.DEG_RIM.to_numpy()
lay = df.layered.to_numpy()
res = {}

# 1. DEG ceiling from containment (big older, small younger)
c = p[(p.kind == "contain") & ~p.near_concentric & ~p.possible_ring]
db, ds = deg[c.big], deg[c.small]
ok = np.isfinite(db) & np.isfinite(ds) & (db != ds)
q = float((ds[ok] > db[ok]).mean())
res["deg_ceiling_from_containment"] = {"q": round(q, 4), "n": int(ok.sum())}


def implied(A):
    return float(np.clip((A - (1 - q)) / (2 * q - 1), 0, 1))


o = p[(p.kind == "overlap") & (p.vote.fillna(0) != 0)].copy()
o["Ds"] = 2 * o.r_small
older = np.where(o.vote > 0, o.big, o.small); younger = np.where(o.vote > 0, o.small, o.big)
do, dy = deg[older], deg[younger]
o["deg_ok"] = np.isfinite(do) & np.isfinite(dy) & (do != dy)
o["deg_agree"] = dy > do
rows = []
for m, g in o.groupby("method"):
    gg = g[g.deg_ok]
    A = float(gg.deg_agree.mean()) if len(gg) else np.nan
    rows.append(dict(method=m, n_edges=len(g), n_with_deg=len(gg), deg_agreement=round(A, 4),
                     implied_accuracy=round(implied(A), 3) if len(gg) else None))
res["deg_check_by_method"] = rows
bins = [1, 3, 5, 10, 20, 1e4]
o["Dbin"] = pd.cut(o.Ds, bins, right=False).astype(str)
rows = []
for (m, b), g in o.groupby(["method", "Dbin"]):
    gg = g[g.deg_ok]
    if len(gg) >= 20:
        A = float(gg.deg_agree.mean())
        rows.append(dict(method=m, D_small=b, n=len(gg), deg_agreement=round(A, 4),
                         implied_accuracy=round(implied(A), 3)))
res["deg_check_by_size"] = rows

# 2. topo vs arc cross agreement
ov = p[p.kind == "overlap"].copy(); ov["Ds"] = 2 * ov.r_small
both = ov[(ov.topo_vote != 0) & (ov.arc_vote != 0)]
res["topo_vs_arc"] = {"n": int(len(both)),
                      "agreement": round(float((both.topo_vote == both.arc_vote).mean()), 4)}

# 3. layered ejecta
lo, ly = lay[older], lay[younger]
one = lo ^ ly
res["layered_ejecta_check"] = {"n": int(one.sum()),
                               "younger_is_the_layered_one": round(float(ly[one].mean()), 4)}

# 4. CTX visual
if a.ctx and os.path.exists(a.ctx):
    v = pd.read_csv(a.ctx)
    v = v.drop(columns=[c for c in ("vote", "method") if c in v.columns])
    cur = p[["vote", "method", "kind"]].rename(columns={"kind": "kind_now"})
    v = v.merge(cur, left_on="pair_id", right_index=True, how="left")
    v.loc[v.kind_now == "contain", "vote"] = 1
    v.loc[v.kind_now == "contain", "method"] = "containment"
    d = v[v.human.notna() & (v.vote.fillna(0) != 0)]
    res["ctx_visual"] = {
        "pairs_shown": int(len(v)),
        "judgeable": int(v.human.notna().sum()),
        "oriented_and_judgeable": int(len(d)),
        "agreement": round(float((d.human == d.vote).mean()), 4),
        "by_method": {m: {"n": int(len(g)), "agree": int((g.human == g.vote).sum())}
                      for m, g in d.groupby("method")},
    }
json.dump(res, open(os.path.join(a.out, "validation.json"), "w"), indent=2)
print(json.dumps(res, indent=2))
