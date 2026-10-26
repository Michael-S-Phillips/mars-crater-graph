"""Evidence for which crater of an overlapping pair is younger.

Sign convention for every vote: +1 = the SMALLER crater is younger (the larger
crater is older), -1 = the LARGER crater is younger, 0 = no decision.

E1  rim continuity on topography (topo_*)
    The younger crater's rim stays a ridge where it crosses the older crater;
    the older crater's rim is erased there. For each crater we measure rim
    prominence p(az) = crest - mean(inner wall at 0.7 r, outer flank at 1.35 r)
    along every azimuth, then compare the median prominence on the arc that lies
    inside the partner ("intrusion arc") with the median on arcs well outside the
    partner ("reference arc"). Preservation P = p_in / p_ref is ~1 for an intact
    (younger) rim and ~0 for a breached (older) rim.

E2  rim-arc completeness from the Robbins catalog (arc_*)
    ARC_IMG is the fraction of the rim Robbins could trace in THEMIS daytime IR.
    If crater k is the older one, the part of its rim inside the partner (fraction
    f_k) should be missing, so 1 - ARC_k >= ~f_k. Youth-likeness
    y_k = 1 - clip((1 - ARC_k) / f_k, 0, 1) is ~1 for an intact rim, ~0 for a
    rim missing exactly where the partner sits. Only defined when f_k >= 0.1.
    A "larger crater is younger" call requires y defined for both craters.

E3  degradation state (deg_*), DEG_RIM 1 (degraded) .. 4 (fresh). Kept separate and
    used for validation, not for orienting edges.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .catalog import R_MARS_KM, destination, bearing

RHO = np.array([0.70, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.35])
CREST = slice(1, 8)


def rim_profiles(df: pd.DataFrame, nodes, dem, px_km: float = 0.463,
                 max_az: int = 720, min_az: int = 72):
    """Rim prominence around each crater in `nodes`. Returns dict node -> arrays."""
    out = {}
    lat = df["lat"].to_numpy(); lon = df["lon"].to_numpy(); r = df["r_km"].to_numpy()
    for n in nodes:
        if abs(lat[n]) > dem.lat_max - 1.0:
            continue
        n_az = int(np.clip(2 * np.pi * r[n] / (0.5 * px_km), min_az, max_az))
        az = np.linspace(0, 2 * np.pi, n_az, endpoint=False)
        ang = (RHO[None, :] * r[n] / R_MARS_KM) * np.ones((n_az, 1))
        la, lo = destination(lat[n], lon[n], az[:, None], ang)
        h = dem.sample(la, lo)
        crest = h[:, CREST].max(axis=1)
        prom = crest - 0.5 * (h[:, 0] + h[:, -1])
        out[n] = dict(az=az, lat=la[:, 4], lon=lo[:, 4], prom=prom.astype(np.float32))
    return out


def _arc_stats(prof, cen_lat, cen_lon, r_partner):
    """Median prominence inside / well outside the partner circle."""
    # angular distance from rim points to the partner center
    la1, lo1 = np.radians(prof["lat"]), np.radians(prof["lon"])
    la2, lo2 = np.radians(cen_lat), np.radians(cen_lon)
    cosang = (np.sin(la1) * np.sin(la2)
              + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2))
    dist = np.arccos(np.clip(cosang, -1, 1)) * R_MARS_KM
    inside = dist < 0.85 * r_partner
    outside = dist > 1.25 * r_partner
    p = prof["prom"]
    p_in = float(np.median(p[inside])) if inside.sum() >= 3 else np.nan
    p_ref = float(np.median(p[outside])) if outside.sum() >= 3 else np.nan
    return p_in, p_ref, int(inside.sum()), int(outside.sum())


def topo_test(pairs: pd.DataFrame, df: pd.DataFrame, profiles: dict,
              min_ref_m: float = 15.0, intact: float = 0.6,
              breached: float = 0.25, min_small_D_km: float = 3.0) -> pd.DataFrame:
    """Two-sided rim-continuity vote. Both rims must have a measurable reference
    prominence (>= min_ref_m) and the smaller crater must be resolved by the DEM
    (D >= min_small_D_km; ~6.5 MOLA 128 ppd pixels). One-sided calls (only one
    rim measurable) are recorded in `topo_single` for diagnostics but are NOT
    used: an unresolved partner cannot visibly breach a rim, so they are biased
    toward calling the resolved crater younger (they disagreed with E2 ~70 % of
    the time in testing)."""
    lat = df["lat"].to_numpy(); lon = df["lon"].to_numpy(); r = df["r_km"].to_numpy()
    rows = []
    for b, s in zip(pairs["big"].to_numpy(), pairs["small"].to_numpy()):
        rec = {}
        for tag, k, m in (("big", b, s), ("small", s, b)):
            prof = profiles.get(k)
            if prof is None:
                rec.update({f"pin_{tag}": np.nan, f"pref_{tag}": np.nan,
                            f"nin_{tag}": 0, f"P_{tag}": np.nan})
                continue
            p_in, p_ref, n_in, n_out = _arc_stats(prof, lat[m], lon[m], r[m])
            P = p_in / p_ref if (np.isfinite(p_in) and np.isfinite(p_ref)
                                 and p_ref >= min_ref_m) else np.nan
            rec.update({f"pin_{tag}": p_in, f"pref_{tag}": p_ref,
                        f"nin_{tag}": n_in, f"P_{tag}": P})
        rows.append(rec)
    t = pd.DataFrame(rows, index=pairs.index)

    def state(P):
        return np.where(~np.isfinite(P), 0, np.where(P >= intact, 1,
                        np.where(P <= breached, -1, 0)))
    sb, ss = state(t.P_big.to_numpy()), state(t.P_small.to_numpy())
    # +1: small intact & big breached; -1: big intact & small breached
    resolved = (2 * pairs["r_small"].to_numpy()) >= min_small_D_km
    vote = np.zeros(len(t), int)
    vote[(ss == 1) & (sb == -1) & resolved] = 1
    vote[(sb == 1) & (ss == -1) & resolved] = -1
    single = np.zeros(len(t), int)  # one side undefined, other decisive
    nb, ns = ~np.isfinite(t.P_big.to_numpy()), ~np.isfinite(t.P_small.to_numpy())
    single[nb & (ss == 1)] = 1; single[nb & (ss == -1)] = -1
    single[ns & (sb == 1)] = -1; single[ns & (sb == -1)] = 1
    t["topo_vote"] = vote
    t["topo_single"] = single
    t["topo_delta"] = t.P_small - t.P_big
    t["topo_conflict"] = resolved & (((ss == 1) & (sb == 1)) | ((ss == -1) & (sb == -1)))
    return t


def arc_test(pairs: pd.DataFrame, df: pd.DataFrame, min_f: float = 0.10,
             thresh: float = 0.5) -> pd.DataFrame:
    arc = df["ARC_IMG"].to_numpy()
    miss_b = 1 - arc[pairs["big"].to_numpy()]
    miss_s = 1 - arc[pairs["small"].to_numpy()]
    fb = pairs["f_big_in_small"].to_numpy(); fs = pairs["f_small_in_big"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        yb = np.where(fb >= min_f, 1 - np.clip(miss_b / fb, 0, 1), np.nan)
        ys = np.where(fs >= min_f, 1 - np.clip(miss_s / fs, 0, 1), np.nan)
    # compare youth-likeness; a missing side is treated as neutral (0.5)
    yb_ = np.where(np.isfinite(yb), yb, 0.5); ys_ = np.where(np.isfinite(ys), ys, 0.5)
    score = ys_ - yb_
    score[~np.isfinite(yb) & ~np.isfinite(ys)] = np.nan
    vote = np.where(score >= thresh, 1, np.where(score <= -thresh, -1, 0))
    vote[~np.isfinite(score)] = 0
    # "larger crater is younger" needs two-sided evidence: the larger rim must be
    # measurably complete where it crosses the smaller one. One-sided -1 calls
    # (smaller rim incomplete, larger rim barely overlapping) agreed with
    # degradation states only ~70 % of the time and failed CTX spot checks.
    vote[(vote == -1) & ~np.isfinite(yb)] = 0
    return pd.DataFrame({"y_big": yb, "y_small": ys, "arc_score": score,
                         "arc_vote": vote}, index=pairs.index)


def deg_compare(pairs: pd.DataFrame, df: pd.DataFrame) -> pd.Series:
    """+1 if the smaller crater has the fresher rim, -1 if the larger does."""
    deg = df["DEG_RIM"].to_numpy()
    db = deg[pairs["big"].to_numpy()]; ds = deg[pairs["small"].to_numpy()]
    v = np.sign(ds - db)
    v[~(np.isfinite(db) & np.isfinite(ds))] = np.nan
    return pd.Series(v, index=pairs.index, name="deg_vote")


def combine(t: pd.DataFrame) -> pd.DataFrame:
    """Merge E1 (two-sided topography) and E2 (rim-arc completeness).

    high      E1 and E2 agree
    medium    E2 decisive, E1 silent or not applicable
    low       E1 decisive, E2 silent (excluded from the default graph)
    conflict  E1 and E2 disagree -> left unoriented
    none      neither is decisive
    """
    tv, av = t.topo_vote.to_numpy(), t.arc_vote.to_numpy()
    vote = np.zeros(len(t), int)
    conf = np.full(len(t), "none", dtype=object)
    method = np.full(len(t), "", dtype=object)
    m = (tv != 0) & (av == tv)
    vote[m] = tv[m]; conf[m] = "high"; method[m] = "topo+arc"
    # topography alone is kept but marked low: where the rim-arc test is silent,
    # the topographic call agreed with degradation states far less often
    # (implied accuracy ~0.8) than when both tests speak (~1.0)
    m = (tv != 0) & (av == 0)
    vote[m] = tv[m]; conf[m] = "low"; method[m] = "topo"
    m = (tv == 0) & (av != 0)
    vote[m] = av[m]; conf[m] = "medium"; method[m] = "arc"
    m = (tv != 0) & (av != 0) & (av != tv)
    vote[m] = 0; conf[m] = "conflict"; method[m] = "conflict"
    return pd.DataFrame({"vote": vote, "confidence": conf, "method": method}, index=t.index)
