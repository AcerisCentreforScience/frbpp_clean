"""Scan a newly released batch on the G4b grid and classify pairs.

Usage: python scripts/scan_newbatch_candidates.py
Output: data/processed/upgrade/newbatch_pairs_<YYYYMMDD>.json
"""
import datetime as dt
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import (EventTable, build_event_table,  # noqa: E402
                                      dm_shuffle_counts, geometry_pairs,
                                      same_event_mask)

THETAS = (60.0, 180.0, 600.0)
DM_WINDOWS = (1.0, 5.0, 20.0)
N_SHUFFLE = 200
rng = np.random.default_rng(20261008)

c1 = pd.read_parquet(ROOT / "data" / "interim" / "cat1_bursts.parquet")
c2 = pd.read_parquet(ROOT / "data" / "interim" / "cat2_bursts.parquet")
ev_cat = build_event_table(c1, c2)

gcn_files = sorted((ROOT / "data" / "raw" / "gcn").glob("history_*.jsonl"))
gcn_path = gcn_files[-1]
vals = []
for ln in gcn_path.read_text(encoding="utf-8").splitlines():
    v = json.loads(json.loads(ln)["value"])
    if v.get("alert_type") == "retraction" or v.get("ra") is None:
        continue
    vals.append(v)
seen, gcn = set(), []
for v in vals:
    if v["id"] not in seen:
        seen.add(v["id"])
        gcn.append(v)

ap = json.loads((ROOT / "data" / "raw" / "apertif" / "apertif_frbs.json").read_text(encoding="utf-8"))
ap_rows = ap["rows"]
n_new = len(gcn) + len(ap_rows)

ra_n = np.concatenate([np.array([v["ra"] for v in gcn]),
                       np.array([r["RAJ2000"] for r in ap_rows])])
dec_n = np.concatenate([np.array([v["dec"] for v in gcn]),
                        np.array([r["DEJ2000"] for r in ap_rows])])
mjd_n = np.concatenate([np.atleast_1d(np.asarray(Time([v["trigger_time"] for v in gcn],
                                                     scale="utc").mjd, dtype=float)),
                        np.array([r["MJD"] for r in ap_rows])])
dm_n = np.concatenate([np.array([v["dm"] for v in gcn], dtype=float),
                       np.array([r["DM"] for r in ap_rows], dtype=float)])
name_n = np.array([f"GCN-{v['id']}" for v in gcn] + [str(r["TNS"]) for r in ap_rows],
                  dtype=object)
cat_n = np.array(["gcn"] * len(gcn) + ["apertif"] * len(ap_rows), dtype=object)
rep_n = np.array([str(v.get("known_source") or "") for v in gcn] + [""] * len(ap_rows),
                 dtype=object)

ev = EventTable(*[np.concatenate([a, b]) for a, b in zip(
    (ev_cat.ra, ev_cat.dec, ev_cat.mjd, ev_cat.dm, ev_cat.name, ev_cat.catalog,
     ev_cat.src_rep, ev_cat.multi_component, ev_cat.excluded),
    (ra_n, dec_n, mjd_n, dm_n, name_n, cat_n, rep_n,
     np.zeros(n_new, bool), np.zeros(n_new, bool)))])
print(f"数据: {gcn_path.name}（GCN {len(gcn)} 条非撤回）+ Apertif {len(ap_rows)}；"
      f"合并表 {len(ev)}（旧 {len(ev_cat)} + 新 {n_new}）\n")

is_newside = np.array([c in ("gcn", "apertif") for c in ev.catalog])


def sep_arcsec(a, b):
    r1, d1, r2, d2 = map(np.radians, [ev.ra[a], ev.dec[a], ev.ra[b], ev.dec[b]])
    x = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    return float(np.degrees(np.arccos(x)) * 3600.0)


rows = []
print("=== 网格：新批次参与配对（分类计数 vs 各自背景）===")
for th in THETAS:
    gi, gj = geometry_pairs(ev, th)
    keep = ~same_event_mask(ev, gi, gj)
    gi, gj = gi[keep], gj[keep]
    m_new = is_newside[gi] | is_newside[gj]
    gi, gj = gi[m_new], gj[m_new]
    for dw in DM_WINDOWS:
        ddm = np.abs(ev.dm[gi] - ev.dm[gj])
        ok = np.isfinite(ddm) & (ddm < dw)
        ia, ib = gi[ok], gj[ok]
        bucket = []
        for a, b in zip(ia.tolist(), ib.tolist()):
            if not is_newside[a]:
                a, b = b, a
            known = ev.src_rep[a] if is_newside[a] else ""
            repb = ev.src_rep[b] if not is_newside[b] else ""
            if ev.catalog[a] in ("gcn", "apertif") and known and not is_newside[b]:
                bucket.append(("A_renewal", a, b))
            elif ev.catalog[a] == "gcn" and not known and repb:
                bucket.append(("B_unassoc_rep", a, b))
            else:
                bucket.append(("C_oneoff", a, b))
        counts = defaultdict(int)
        for tag, _a, _b in bucket:
            counts[tag] += 1
        bgline = []
        for tag in ("A_renewal", "B_unassoc_rep", "C_oneoff"):
            idx = [k for k, (t, _a, _b) in enumerate(bucket) if t == tag]
            if not idx:
                bgline.append(f"{tag}=0")
                continue
            sub_i = np.array([ia[k] for k in idx])
            sub_j = np.array([ib[k] for k in idx])
            bg = dm_shuffle_counts(ev, sub_i, sub_j, dw, 100, rng)
            mean = float(bg.mean())
            thr = mean + 3.0 * np.sqrt(mean)
            flag = "***" if len(idx) > thr else ""
            bgline.append(f"{tag}={len(idx)}(bg {mean:.2f}{flag})")
        print(f"θ={th:4.0f}″ DM±{dw:4.0f}: " + " | ".join(bgline))
        for tag, a, b in bucket:
            rows.append({"theta": th, "dmw": dw, "tag": tag,
                         "new": str(ev.name[a]), "newcat": str(ev.catalog[a]),
                         "known": str(ev.src_rep[a]), "new_mjd": float(ev.mjd[a]),
                         "old": str(ev.name[b]), "oldcat": str(ev.catalog[b]),
                         "old_rep": str(ev.src_rep[b]), "old_mjd": float(ev.mjd[b]),
                         "sep": sep_arcsec(a, b),
                         "ddm": float(abs(ev.dm[a] - ev.dm[b]))})

print("\n=== A 续暴汇总（新事件 → 历史源）===")
agg = defaultdict(lambda: {"n": 0, "old_names": set(), "sep": [], "ddm": [], "reps": set()})
for r in rows:
    if r["tag"] != "A_renewal":
        continue
    g = agg[(r["new"], r["known"])]
    g["n"] += 1
    g["old_names"].add(r["old"])
    g["sep"].append(r["sep"])
    g["ddm"].append(r["ddm"])
    g["reps"].add(r["old_rep"])
for (new, known), g in sorted(agg.items()):
    print(f"  {new} (known={known}) → 历史行 {len(g['old_names'])} 个、配对 {g['n']} 次 "
          f"| sep {min(g['sep']):.0f}–{max(g['sep']):.0f}″ | ΔDM {min(g['ddm']):.2f}–{max(g['ddm']):.2f} "
          f"| reps={sorted(g['reps'])}")

print("\n=== B 未关联 × 重复源 行（唯一对）===")
seen_b = set()
for r in rows:
    if r["tag"] != "B_unassoc_rep" or (r["new"], r["old"]) in seen_b:
        continue
    seen_b.add((r["new"], r["old"]))
    print(f"  {r['new']} × {r['old']}[rep={r['old_rep']}] sep={r['sep']:.0f}″ "
          f"ΔDM={r['ddm']:.2f} Δt={abs(r['new_mjd']-r['old_mjd']):.0f}d")
if not seen_b:
    print("  无")

print("\n=== C one-off 候选（唯一对；新×旧 or 新×新）===")
seen_c = set()
for r in sorted(rows, key=lambda x: (x["ddm"], x["sep"])):
    if r["tag"] != "C_oneoff" or (r["new"], r["old"]) in seen_c:
        continue
    seen_c.add((r["new"], r["old"]))
    print(f"  {r['new']} ({r['newcat']}) × {r['old']} ({r['oldcat']}) "
          f"sep={r['sep']:.0f}″ ΔDM={r['ddm']:.2f} Δt={abs(r['new_mjd']-r['old_mjd']):.0f}d")
if not seen_c:
    print("  无")

outdir = ROOT / "data" / "processed" / "upgrade"
outdir.mkdir(parents=True, exist_ok=True)
out = outdir / f"newbatch_pairs_{dt.datetime.now():%Y%m%d}.json"
out.write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"\n配对明细 {len(rows)} 条 → {out.relative_to(ROOT)}")