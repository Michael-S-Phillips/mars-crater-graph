"""Extract CTX (Murray Lab global mosaic, 5 m/px, 2x2 degree tiles) chips around
crater-pair windows for visual validation.

  python extract_ctx_chips.py windows.json OUTDIR CTX_ROOT
CTX_ROOT holds folders like W28_S32/CTX-Mosaic_MurrayLab_Feb2018_W26_S30.tif, where a
tile named W<a>_S<b> spans longitude -a..-a+2 and latitude -b..-b+2.
"""
import glob, json, math, os, sys
import numpy as np
import rasterio
from rasterio.merge import merge

R_M = 3396190.0
K = math.radians(1) * R_M
wins = json.load(open(sys.argv[1])); out = sys.argv[2]; root = sys.argv[3]
os.makedirs(out, exist_ok=True)
index = {}
for f in glob.glob(os.path.join(root, "W*", "*.tif")):
    nm = os.path.basename(f).replace(".tif", "").split("_")
    index[(int(nm[-2][1:]), int(nm[-1][1:]))] = f

N = 480
for w in wins:
    f_out = os.path.join(out, f"{w['pair_id']}.npz")
    if os.path.exists(f_out):
        continue
    lon = w["clon"] - 360 if w["clon"] > 180 else w["clon"]
    lat = w["clat"]; hw = w["hw_km"]
    dlat = math.degrees(hw * 1000 / R_M); dlon = dlat / math.cos(math.radians(lat))
    W, S, E, Nn = lon - dlon, lat - dlat, lon + dlon, lat + dlat
    tiles = set()
    for lo in (W, E):
        for la in (S, Nn):
            a = -2 * math.floor(lo / 2); b = -2 * math.floor(la / 2)
            if (a, b) in index:
                tiles.add(index[(a, b)])
    if not tiles:
        continue
    dss = [rasterio.open(t) for t in sorted(tiles)]
    res = (2 * dlon * K / N, 2 * dlat * K / N)
    arr, _ = merge(dss, bounds=(W * K, S * K, E * K, Nn * K), res=res, nodata=0)
    for d in dss:
        d.close()
    np.savez_compressed(f_out, ctx=arr[0, :N, :N], **{k: np.float64(v) for k, v in w.items()})
print("ok")
