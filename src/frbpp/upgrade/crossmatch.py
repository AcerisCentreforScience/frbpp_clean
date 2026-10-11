"""Event-table construction, per-burst spatio-temporal cones, DM-window clustering, cross-catalogue naming."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

SAME_EVENT_DT_DAYS = 0.05
LOGIC_VERSION = "1.3.0"
LOGIC_VERSION_V12 = "1.2.0"
_TRUEISH = {"True", "true", "1", "1.0"}


@dataclass
class EventTable:
    ra: np.ndarray
    dec: np.ndarray
    mjd: np.ndarray
    dm: np.ndarray
    name: np.ndarray
    catalog: np.ndarray
    src_rep: np.ndarray
    multi_component: np.ndarray
    excluded: np.ndarray

    def __len__(self):
        return int(self.ra.size)


def _num(s) -> np.ndarray:
    return pd.to_numeric(s, errors="coerce").to_numpy(dtype=float)


def _str(s) -> np.ndarray:
    return s.astype(str).str.strip().to_numpy()


def clean_rep(s) -> np.ndarray:
    rep = s.astype("string").str.strip()
    rep = rep.fillna("")
    rep = rep.mask(rep.isin(["-9999", "nan", "NaN", "None", "none", ""]), "")
    return rep.to_numpy(dtype=object)


def build_event_table(cat1_df: pd.DataFrame, cat2_df: pd.DataFrame,
                      event_dedup: bool = True) -> EventTable:
    rep1 = clean_rep(cat1_df.repeater_name)
    rep2 = clean_rep(cat2_df.repeater_name)
    ra = np.concatenate([_num(cat1_df.ra), _num(cat2_df.ra)])
    dec = np.concatenate([_num(cat1_df.dec), _num(cat2_df.dec)])
    mjd = np.concatenate([_num(cat1_df.mjd_400), _num(cat2_df.mjd_400)])
    dm = np.concatenate([_num(cat1_df.dm_fitb), _num(cat2_df.dm_fitb)])
    name = np.concatenate([_str(cat1_df.tns_name), _str(cat2_df.tns_name)])
    cat = np.concatenate([np.full(len(cat1_df), "cat1"), np.full(len(cat2_df), "cat2")])
    rep = np.concatenate([rep1, rep2])
    multi = np.concatenate([_num(cat1_df.sub_num) >= 1,
                            np.zeros(len(cat2_df), dtype=bool)])
    ex1 = cat1_df.excluded_flag.astype(str).isin(_TRUEISH).to_numpy()
    ex2 = (cat2_df.excluded_flag.astype(str).isin(_TRUEISH)
           | cat2_df.sidelobe_flag.astype(str).isin(_TRUEISH)).to_numpy()
    ex = np.concatenate([ex1, ex2])
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(mjd)
    if event_dedup:
        ok = ok & _event_layer_keep(name, cat, mjd, ok)
    return EventTable(ra[ok], dec[ok], mjd[ok], dm[ok], name[ok], cat[ok],
                      rep[ok], multi[ok], ex[ok])


def _event_layer_keep(name: np.ndarray, cat: np.ndarray, mjd: np.ndarray,
                      usable: np.ndarray) -> np.ndarray:
    from collections import defaultdict
    c2mj: dict[str, list] = defaultdict(list)
    for k in range(name.size):
        if cat[k] == "cat2" and usable[k]:
            c2mj[name[k]].append(mjd[k])
    keep = np.ones(name.size, dtype=bool)
    for k in range(name.size):
        if cat[k] == "cat1":
            for t in c2mj.get(name[k], ()):
                if abs(mjd[k] - t) < SAME_EVENT_DT_DAYS:
                    keep[k] = False
                    break
    return keep


def geometry_pairs(ev: EventTable, theta_arcsec: float):
    coord = SkyCoord(ev.ra * u.deg, ev.dec * u.deg)
    idx1, idx2, _, _ = search_around_sky(coord, coord, theta_arcsec * u.arcsec)
    keep = idx1 < idx2
    return idx1[keep], idx2[keep]


def same_event_mask(ev: EventTable, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    dt = np.abs(ev.mjd[i] - ev.mjd[j])
    cross = ev.catalog[i] != ev.catalog[j]
    rule1 = (ev.name[i] == ev.name[j]) & (dt < SAME_EVENT_DT_DAYS)
    rule2 = cross & (dt < 0.01)
    return rule1 | rule2


def cone_pairs(ev: EventTable, theta_arcsec: float, dm_window: float):
    i, j = geometry_pairs(ev, theta_arcsec)
    ddm = np.abs(ev.dm[i] - ev.dm[j])
    dm_ok = np.isfinite(ddm) & (ddm < dm_window)
    m = (~same_event_mask(ev, i, j)) & dm_ok
    return i[m], j[m]


def dm_shuffle_counts(ev: EventTable, i: np.ndarray, j: np.ndarray,
                      dm_window: float, n_shuffle: int, rng) -> np.ndarray:
    finite = np.isfinite(ev.dm)
    base = ev.dm.copy()
    out = np.empty(n_shuffle, dtype=float)
    for k in range(n_shuffle):
        d = base.copy()
        d[finite] = rng.permutation(d[finite])
        ddm = np.abs(d[i] - d[j])
        out[k] = float(np.sum(np.isfinite(ddm) & (ddm < dm_window)))
    return out


def components(n_events: int, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    parent = np.arange(n_events)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(i.tolist(), j.tolist()):
        ra_, rb_ = find(a), find(b)
        if ra_ != rb_:
            parent[ra_] = rb_
    return np.array([find(x) for x in range(n_events)])


def candidate_groups(ev: EventTable, i, j, min_events: int = 2):
    if i.size == 0:
        return []
    roots = components(len(ev), i, j)
    touched = np.unique(np.concatenate([i, j]))
    out = []
    for r in np.unique(roots[touched]):
        members = touched[roots[touched] == r]
        if members.size >= min_events:
            out.append(members)
    return out


def repeater_recall(ev: EventTable, i, j):
    got = {}
    for a, b in zip(i.tolist(), j.tolist()):
        ra_, rb_ = ev.src_rep[a], ev.src_rep[b]
        if ra_ and ra_ == rb_:
            got[ra_] = got.get(ra_, 0) + 1
    known = sorted({r for r in ev.src_rep if r})
    recovered = [r for r in known if got.get(r, 0) >= 1]
    return known, recovered


def name_census(cat1_df: pd.DataFrame, cat2_df: pd.DataFrame) -> dict:
    p1 = cat1_df[cat1_df.is_burst_primary].copy()
    p1["rep1"] = clean_rep(p1.repeater_name)
    c2 = cat2_df.copy()
    c2["mjd_n"] = _num(c2.mjd_400)
    p1["mjd_n"] = _num(p1.mjd_400)
    oneoff = set(p1.loc[p1.rep1 == "", "tns_name"].astype(str))
    last1 = p1.groupby("tns_name")["mjd_n"].max()
    upgrades = []
    for n in sorted(oneoff & set(c2.tns_name.astype(str))):
        tmax1 = float(last1[n])
        new = c2[(c2.tns_name.astype(str) == n) & (c2.mjd_n > tmax1 + 1.0)]
        if len(new):
            upgrades.append({"name": n, "cat1_last_mjd": tmax1,
                             "cat2_new_mjd": float(new.mjd_n.min()),
                             "n_new": int(len(new))})
    return {"n_cat1_oneoff": len(oneoff), "upgrades": upgrades}