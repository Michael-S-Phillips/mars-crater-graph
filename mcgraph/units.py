"""Host geologic unit for each crater from the Tanaka et al. (2014) global map
(USGS SIM 3292). Used as context (which terrain a crater sits on), not as
stratigraphic edges: a crater's host unit only bounds its age when the crater is
superposed on that unit rather than embayed by it, which the map cannot tell us.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PERIOD_ORDER = {"N": 0, "H": 1, "A": 2}
EPOCH_ORDER = {"eN": 0, "mN": 1, "lN": 2, "eH": 3, "lH": 4, "eA": 5, "mA": 6, "lA": 7}


def parse_unit(code: str):
    """Return (oldest period letter, oldest epoch or None) from a SIM 3292 unit code."""
    if not isinstance(code, str) or not code:
        return None, None
    if code[0] in "eml" and len(code) > 1 and code[1] in "NHA":
        return code[1], code[:2]
    if code[:2] in ("AH", "HN", "AN"):
        return code[1], None
    if code[0] in "NHA":
        return code[0], None
    return None, None


def host_units(df: pd.DataFrame, shp_path: str) -> pd.DataFrame:
    import geopandas as gpd
    g = gpd.read_file(shp_path)[["Unit", "UnitDesc", "geometry"]]
    lon = np.where(df.lon > 180, df.lon - 360, df.lon)
    pts = gpd.GeoDataFrame({"node": df.index.to_numpy()},
                           geometry=gpd.points_from_xy(lon, df.lat), crs=g.crs)
    j = gpd.sjoin(pts, g, how="left", predicate="within")
    j = j.drop_duplicates("node").set_index("node")
    out = pd.DataFrame(index=df.index)
    out["unit"] = j["Unit"]
    out["unit_desc"] = j["UnitDesc"]
    pe = out["unit"].map(lambda u: parse_unit(u))
    out["unit_period"] = pe.map(lambda t: t[0])
    out["unit_epoch"] = pe.map(lambda t: t[1])
    return out
