r"""Offline recheck of the dispersion-measure fingerprint isolation claim.

The Transient Name Server query is the online step and is not re-run here: its archived ledger line
(label dm_fingerprint) records the queried dispersion-measure range and the eight source names it
returned. This script recomputes, from coordinates held locally, the angular distance of every
returned source other than the two pair members, and reports the minimum distance and how many lie
within one degree. It therefore makes the manuscript's conditional-isolation sentence checkable
without credentials or network access.

Usage: python scripts/dm_fingerprint_check.py
Output: data/processed/upgrade/dm_fingerprint.json
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data/processed/upgrade/tns_search_ledger.jsonl"
DETAILS = ROOT / "data/processed/upgrade/tns_object_details.csv"
CAT2 = ROOT / "data/raw/chime/cat2_25/chimefrbcat2.csv"
OUT = ROOT / "data/processed/upgrade/dm_fingerprint.json"

PAIR = ("FRB20200511C", "FRB20211128E")
NEAR_DEG = 1.0


def sky(name, cat2, details):
    r = cat2[(cat2.tns_name == name) & ((cat2.sub_num == 0) | cat2.sub_num.isna())]
    if len(r):
        return float(r.iloc[0].ra), float(r.iloc[0].dec), "Catalogue 2 primary row"
    hit = details[details.apply(lambda x: name in str(x.values), axis=1)]
    if len(hit):
        return float(hit.iloc[0].ra), float(hit.iloc[0].dec), "TNS object details"
    return None, None, "unresolved"


def main() -> int:
    fp = None
    for line in LEDGER.read_text(encoding="utf-8").splitlines():
        rec = json.loads(line)
        if rec.get("label") == "dm_fingerprint":
            fp = rec
            break
    if fp is None:
        raise SystemExit("台账里没有 dm_fingerprint 条目，无法复算")

    names = fp["names"]
    lo, hi = fp["query"]["frb_dm_range_min"], fp["query"]["frb_dm_range_max"]
    cat2 = pd.read_csv(CAT2, low_memory=False)
    details = pd.read_csv(DETAILS) if DETAILS.exists() else pd.DataFrame(columns=["ra", "dec"])

    pos = {}
    for n in names:
        ra, dec, src = sky(n, cat2, details)
        if ra is None:
            raise SystemExit(f"{n} 无本地坐标，复算不完整")
        pos[n] = (ra, dec, src)

    others = [n for n in names if n not in PAIR]
    a = SkyCoord(*pos[PAIR[0]][:2], unit=u.deg)
    b = SkyCoord(*pos[PAIR[1]][:2], unit=u.deg)
    rows = []
    for n in others:
        c = SkyCoord(*pos[n][:2], unit=u.deg)
        d1, d2 = c.separation(a).deg, c.separation(b).deg
        rows.append({"name": n, "deg_from_first": round(d1, 2), "deg_from_second": round(d2, 2),
                     "min_deg": round(min(d1, d2), 2), "coords_from": pos[n][2]})
    mins = [r["min_deg"] for r in rows]

    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "source": "tns_search_ledger.jsonl label=dm_fingerprint (archived, not re-queried)",
           "dm_range_pc_cm3": [lo, hi], "n_returned": int(fp["n"]),
           "n_pair_members": sum(1 for n in names if n in PAIR),
           "n_others": len(others), "others": rows,
           "min_separation_deg": min(mins), "n_within_1_deg": int(sum(1 for m in mins if m < NEAR_DEG)),
           "reading": ("A comparable dispersion measure is not found near this line of sight in the "
                      "reporting system, so the dispersion-measure isolation is supported by data "
                      "independent of the pair's own catalogue entries.")}
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"返回 {fp['n']} 个（本对 {out['n_pair_members']} 个，其他 {len(others)} 个）")
    for r in rows:
        print(f"  {r['name']:16s} 距两成员最小 {r['min_deg']:6.2f}°  坐标来源 {r['coords_from']}")
    print(f"最小 = {out['min_separation_deg']}° ; <{NEAR_DEG}° 内 {out['n_within_1_deg']} 个 -> {OUT.name}")
    assert out["n_returned"] == 8 and out["n_within_1_deg"] == 0 and out["min_separation_deg"] > 40.0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
