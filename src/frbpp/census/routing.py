"""Stratified source routing rules and hygiene-field assignment."""
from __future__ import annotations

import numpy as np
import pandas as pd

from frbpp.pointprocess.decluster import decluster
from frbpp.upgrade.crossmatch import clean_rep

MIN_EPISODES_A = 5
T_CLUSTER_D = 1.0
MJD_BUCKET_D = 1e-6


def build_event_frame(cat2_df: pd.DataFrame) -> pd.DataFrame:
    df = cat2_df.copy()
    rep = clean_rep(df.repeater_name)
    src = np.where(np.array([r != "" for r in rep]), rep, df.tns_name.astype(str).to_numpy())
    df["src"] = src
    df["mjd_n"] = pd.to_numeric(df.mjd_400, errors="coerce").to_numpy(dtype=float)
    df["ra_n"] = pd.to_numeric(df.ra, errors="coerce").to_numpy(dtype=float)
    df["dec_n"] = pd.to_numeric(df.dec, errors="coerce").to_numpy(dtype=float)
    df["dm_n"] = pd.to_numeric(df.dm_fitb, errors="coerce").to_numpy(dtype=float)
    df = df[np.isfinite(df.mjd_n)].copy()
    df["_pri"] = pd.to_numeric(df.sub_num, errors="coerce").fillna(0).eq(0)
    df["_bucket"] = (df.mjd_n / MJD_BUCKET_D).round().astype("int64")
    df = (df.sort_values(["src", "_bucket", "_pri"], ascending=[True, True, False])
            .drop_duplicates(["src", "_bucket"], keep="first"))
    return df[["src", "tns_name", "mjd_n", "ra_n", "dec_n", "dm_n"]].reset_index(drop=True)


def route_sources(events: pd.DataFrame):
    rows, episodes = [], {}
    for src, g in events.groupby("src", sort=True):
        ep, _ = decluster(g.mjd_n.to_numpy(), t_cluster_days=T_CLUSTER_D)
        episodes[src] = ep
        n_ep = int(ep.size)
        has_pos = bool(np.isfinite([g.ra_n.iloc[0], g.dec_n.iloc[0]]).all())
        cat = "A" if n_ep >= MIN_EPISODES_A else ("B" if n_ep >= 2 else ("C" if has_pos else "D"))
        mjd = g.mjd_n.to_numpy()
        span = float(mjd.max() - mjd.min()) if len(mjd) > 1 else 0.0
        mid = float(mjd.min()) + span / 2.0
        r1, r2 = int((mjd < mid).sum()), int((mjd >= mid).sum())
        t1, t2 = max(mid - float(mjd.min()), 1e-9), max(float(mjd.max()) - mid, 1e-9)
        gaps = np.diff(ep) if ep.size > 1 else np.array([])
        rel_scatter = (float(np.median(np.abs(gaps - np.median(gaps))) / np.median(gaps))
                       if gaps.size else np.nan)
        rows.append({
            "src": src, "category": cat, "n_bursts": int(len(g)), "n_episodes": n_ep,
            "mjd_min": float(mjd.min()), "mjd_max": float(mjd.max()), "span_days": span,
            "ra": float(g.ra_n.iloc[0]) if has_pos else np.nan,
            "dec": float(g.dec_n.iloc[0]) if has_pos else np.nan,
            "rate1_per_day": r1 / t1, "rate2_per_day": r2 / t2,
            "relative_scatter": rel_scatter,
            "sigma_inflation": np.nan,
        })
    return pd.DataFrame(rows), episodes