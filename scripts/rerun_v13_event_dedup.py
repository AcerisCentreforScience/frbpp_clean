"""Rerun the G4b grid, candidate union and isolation MC under the v1.3 event-layer de-duplication.

Usage: python scripts/rerun_v13_event_dedup.py
Output: data/processed/upgrade/g4b_rerun_20261009_v13.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import (                       # noqa: E402
    LOGIC_VERSION, LOGIC_VERSION_V12, EventTable, build_event_table,
    candidate_groups, cone_pairs, dm_shuffle_counts, geometry_pairs,
    same_event_mask)

RAW1 = ROOT / "data" / "interim" / "cat1_bursts.parquet"
RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
UP = ROOT / "data" / "processed" / "upgrade"

THETAS = (60.0, 180.0, 600.0)
DMWS = (1.0, 5.0, 20.0)
N_SHUFFLE = 200
SEED = 20261009


def subset(ev: EventTable, m: np.ndarray) -> EventTable:
    return EventTable(*(a[m] for a in (ev.ra, ev.dec, ev.mjd, ev.dm, ev.name,
                                       ev.catalog, ev.src_rep, ev.multi_component,
                                       ev.excluded)))


def grid_run(ev: EventTable, rng, dedup_bg: bool = True):
    rows = []
    for th in THETAS:
        gi, gj = geometry_pairs(ev, th)
        keep = ~same_event_mask(ev, gi, gj)
        for dw in DMWS:
            ddm = np.abs(ev.dm[gi] - ev.dm[gj])
            obs = int(np.sum(keep & np.isfinite(ddm) & (ddm < dw)))
            bgi, bgj = (gi[keep], gj[keep]) if dedup_bg else (gi, gj)
            bg = dm_shuffle_counts(ev, bgi, bgj, dw, N_SHUFFLE, rng)
            mean = float(bg.mean())
            thr = mean + 3.0 * np.sqrt(mean)
            i, j = gi[keep & np.isfinite(ddm) & (ddm < dw)], gj[keep & np.isfinite(ddm) & (ddm < dw)]
            rows.append({"theta": th, "dm_window": dw, "obs": obs, "bg_mean": mean,
                         "thr": thr, "sig": bool(obs > thr),
                         "ngroups": len(candidate_groups(ev, i, j)), "i": i, "j": j})
    return rows


def union_keys(cfgs, ev):
    best = {}
    for c in cfgs:
        if not c["sig"]:
            continue
        for a, b in zip(c["i"].tolist(), c["j"].tolist()):
            key = tuple(sorted([str(ev.name[a]), str(ev.name[b])]))
            if key not in best or c["theta"] < best[key]:
                best[key] = c["theta"]
    return best


def main() -> int:
    c1, c2 = pd.read_parquet(RAW1), pd.read_parquet(RAW2)

    res = {"logic_version": LOGIC_VERSION, "logic_version_prev": LOGIC_VERSION_V12,
           "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "change": "事件层去重（v1.3）：先剔除与 Cat2 同名且 |Δt|<0.05 d 的 Cat1 行，再建模",
           "config": {"thetas": list(THETAS), "dm_windows": list(DMWS),
                      "n_shuffle": N_SHUFFLE, "seed": SEED}}

    out = {}
    for tag, dd in (("v12", False), ("v13", True)):
        ev = build_event_table(c1, c2, event_dedup=dd)
        oo = np.array([r == "" for r in ev.src_rep])
        ev_oo = subset(ev, oo)
        fin = oo & np.isfinite(ev.dm)
        rng = np.random.default_rng(SEED)
        cfgs = grid_run(ev_oo, rng)
        union = union_keys(cfgs, ev_oo)
        n_sig = int(sum(c["sig"] for c in cfgs))
        n_pairs = int(sum(c["obs"] for c in cfgs if c["sig"]))
        rec = {
            "n_events_all": int(len(ev)),
            "n_events_oneoff": int(oo.sum()),
            "n_events_oneoff_finite_dm": int(fin.sum()),
            "n_sig_cfg": n_sig,
            "n_candidate_pairs_union": len(union),
            "n_pairs_in_sig_cfgs": n_pairs,
            "sum_groups_sig": int(sum(c["ngroups"] for c in cfgs if c["sig"])),
            "grid": [{k: c[k] for k in ("theta", "dm_window", "obs", "bg_mean",
                                        "thr", "sig", "ngroups")} for c in cfgs],
            "union_keys": sorted("|".join(k) for k in union),
        }
        out[tag] = rec
        print(f"[{tag}] events {rec['n_events_all']} | one-off {rec['n_events_oneoff']} "
              f"(finite DM {rec['n_events_oneoff_finite_dm']}) | sig cfg {n_sig}/9 | "
              f"union {len(union)} pairs")
    res["v12_row_level"] = out["v12"]
    res["v13_event_dedup"] = out["v13"]

    same = out["v12"]["union_keys"] == out["v13"]["union_keys"]
    res["invariance"] = {
        "candidate_union_identical": bool(same),
        "v12_union_n": len(out["v12"]["union_keys"]),
        "v13_union_n": len(out["v13"]["union_keys"]),
        "only_in_v12": sorted(set(out["v12"]["union_keys"]) - set(out["v13"]["union_keys"])),
        "only_in_v13": sorted(set(out["v13"]["union_keys"]) - set(out["v12"]["union_keys"])),
        "interpretation": ("候选并集逐对比对完全一致 → 去重提前到事件层不改变任何候选、"
                           "p 值或结论；仅以事件池为分母的量（子空间规模/零假设/注入目录）缩小。"
                           if same else "候选并集发生变化——必须逐项复核并更新论文。"),
    }
    diff_cells = []
    for a, b in zip(out["v12"]["grid"], out["v13"]["grid"]):
        if a["obs"] != b["obs"] or a["ngroups"] != b["ngroups"]:
            diff_cells.append({"theta": a["theta"], "dm_window": a["dm_window"],
                               "v12_obs": a["obs"], "v13_obs": b["obs"],
                               "v12_ngroups": a["ngroups"], "v13_ngroups": b["ngroups"]})
    res["invariance"]["cells_with_difference"] = diff_cells

    UP.mkdir(parents=True, exist_ok=True)
    (UP / "g4b_rerun_20261009_v13.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(f"\n→ {UP / 'g4b_rerun_20261009_v13.json'}")
    print(f"候选并集一致: {same}（v12={len(out['v12']['union_keys'])} v13={len(out['v13']['union_keys'])}）")
    print(f"有差异的档位: {len(diff_cells)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
