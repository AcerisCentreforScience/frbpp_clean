"""Stage 4b acceptance gate G4b.

Usage: python scripts/gate_g4b.py
Output: data/processed/upgrade/g4b_report.json
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

from frbpp.upgrade.crossmatch import (            # noqa: E402
    EventTable, build_event_table, candidate_groups, clean_rep, cone_pairs,
    dm_shuffle_counts, geometry_pairs, name_census, repeater_recall, same_event_mask)

RAW1 = ROOT / "data" / "interim" / "cat1_bursts.parquet"
RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
OUT = ROOT / "data" / "processed" / "upgrade"

THETAS = (60.0, 180.0, 600.0)
DM_WINDOWS = (1.0, 5.0, 20.0)
N_SHUFFLE = 200
RECALL_MIN = 40


def subset(ev: EventTable, m: np.ndarray) -> EventTable:
    return EventTable(*(a[m] for a in (ev.ra, ev.dec, ev.mjd, ev.dm, ev.name,
                                       ev.catalog, ev.src_rep, ev.multi_component,
                                       ev.excluded)))


def grid_run(ev: EventTable, rng) -> list[dict]:
    rows = []
    for th in THETAS:
        gi, gj = geometry_pairs(ev, th)
        dedup_keep = ~same_event_mask(ev, gi, gj)
        for dw in DM_WINDOWS:
            ddm = np.abs(ev.dm[gi] - ev.dm[gj])
            obs = int(np.sum(dedup_keep & np.isfinite(ddm) & (ddm < dw)))
            bg = dm_shuffle_counts(ev, gi, gj, dw, N_SHUFFLE, rng)
            mean = float(bg.mean())
            thr = mean + 3.0 * np.sqrt(mean)
            i, j = cone_pairs(ev, th, dw)
            rows.append({"theta_arcsec": th, "dm_window": dw, "n_geom": int(gi.size),
                         "observed_pairs": obs, "background_mean": mean,
                         "background_p99": float(np.quantile(bg, 0.99)),
                         "threshold_3sig": thr, "significant": bool(obs > thr),
                         "n_groups": len(candidate_groups(ev, i, j))})
    return rows


def main() -> int:
    rng = np.random.default_rng(20261010)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G4b 冻结口径 v1 (+v1.1 操作化)"}

    c1 = pd.read_parquet(RAW1)
    c2 = pd.read_parquet(RAW2)
    ev = build_event_table(c1, c2)
    print(f"事件表: {len(ev)} 行（剔除无坐标 {len(c1) + len(c2) - len(ev)}）| "
          f"重复源事件 {int(sum(1 for r in ev.src_rep if r))} | "
          f"排除标记 {int(ev.excluded.sum())} | 子暴行 {int(ev.multi_component.sum())}")
    out["events"] = {"n": len(ev), "n_dropped_no_coords": len(c1) + len(c2) - len(ev),
                     "n_rep_events": int(sum(1 for r in ev.src_rep if r)),
                     "n_excluded_flagged": int(ev.excluded.sum()),
                     "n_sub_burst_rows": int(ev.multi_component.sum())}

    cen = name_census(c1, c2)
    f1 = c2[c2.catalog1_flag.astype(str) == "1"]
    suffix = sorted(set(f1.tns_name.astype(str)) - set(c1.tns_name.astype(str)))
    print(f"[census] Cat1 one-off {cen['n_cat1_oneoff']} 源 | 名称基升级 = "
          f"{len(cen['upgrades'])}（期望 0）| 异后缀疑似同源 {suffix}")
    out["name_census"] = cen
    out["suffix_variant_names"] = suffix

    full_rows = grid_run(ev, rng)
    for r in full_rows:
        print(f"  [全域] θ={r['theta_arcsec']:5.0f}\" DM±{r['dm_window']:4.0f}: "
              f"观测 {r['observed_pairs']:5d} | 背景 {r['background_mean']:8.1f} | 组 {r['n_groups']}")
    out["grid_full"] = full_rows

    i6, j6 = cone_pairs(ev, 600.0, 5.0)
    known, recovered = repeater_recall(ev, i6, j6)
    ok1 = len(recovered) >= RECALL_MIN
    print(f"[G4b①] 重复源召回 {len(recovered)}/{len(known)}（判 ≥{RECALL_MIN}）→ "
          f"{'OK' if ok1 else 'FAIL'}")
    out["g4b1"] = {"n_known": len(known), "n_recovered": len(recovered),
                   "recovered": recovered, "criterion": RECALL_MIN, "ok": bool(ok1)}

    ev_oo = subset(ev, np.array([r == "" for r in ev.src_rep]))
    print(f"\none-off 子空间: {len(ev_oo)} 事件")
    oo_rows = grid_run(ev_oo, rng)
    for r in oo_rows:
        print(f"  [one-off] θ={r['theta_arcsec']:5.0f}\" DM±{r['dm_window']:4.0f}: "
              f"观测 {r['observed_pairs']:4d} | 背景 {r['background_mean']:6.2f} | "
              f"显著 {r['significant']} | 组 {r['n_groups']}")
    out["grid_oneoff"] = oo_rows

    best = {}
    for r in oo_rows:
        if not r["significant"]:
            continue
        i, j = cone_pairs(ev_oo, r["theta_arcsec"], r["dm_window"])
        for a, b in zip(i.tolist(), j.tolist()):
            key = tuple(sorted([str(ev_oo.name[a]), str(ev_oo.name[b])]))
            cur = {"theta": r["theta_arcsec"], "dm_window": r["dm_window"],
                   "name1": str(ev_oo.name[a]), "cat1": str(ev_oo.catalog[a]),
                   "name2": str(ev_oo.name[b]), "cat2": str(ev_oo.catalog[b]),
                   "dm1": float(ev_oo.dm[a]), "dm2": float(ev_oo.dm[b]),
                   "dt_days": float(abs(ev_oo.mjd[a] - ev_oo.mjd[b]))}
            if key not in best or cur["theta"] < best[key]["theta"]:
                best[key] = cur
    cand = sorted(best.values(), key=lambda c: (c["dt_days"], c["theta"]))
    same_day = [c for c in cand if c["dt_days"] < 1.0]
    cross_day = [c for c in cand if c["dt_days"] >= 1.0]
    print(f"[候选] 唯一配对 {len(cand)} = 同日档 {len(same_day)} + 跨日档 {len(cross_day)}")
    for c in cross_day[:10]:
        print(f"   [跨日] {c['name1']} × {c['name2']} θ={c['theta']:.0f}\" "
              f"ΔDM={abs(c['dm1'] - c['dm2']):.1f} Δt={c['dt_days']:.1f}d")
    out["candidates"] = {"same_day": same_day, "cross_day": cross_day}

    n_sig_ctrl = 0
    ctrl_pairs = []
    for _ in range(20):
        ra_p = rng.permutation(ev_oo.ra)
        ev_c = EventTable(ra_p, ev_oo.dec, ev_oo.mjd, ev_oo.dm, ev_oo.name,
                          ev_oo.catalog, ev_oo.src_rep, ev_oo.multi_component,
                          ev_oo.excluded)
        ic, jc = cone_pairs(ev_c, 180.0, 5.0)
        obs_c = int(ic.size)
        gi, gj = geometry_pairs(ev_c, 180.0)
        bg_c = dm_shuffle_counts(ev_c, gi, gj, 5.0, 100, rng)
        thr_c = float(bg_c.mean()) + 3.0 * np.sqrt(bg_c.mean())
        ctrl_pairs.append((obs_c, float(bg_c.mean())))
        n_sig_ctrl += int(obs_c > thr_c)
    obs_ctl = np.array([p[0] for p in ctrl_pairs], dtype=float)
    ok2 = bool(n_sig_ctrl <= 1)
    print(f"[G4b②] 逐行 RA 置换 20 次: 对数 均值={obs_ctl.mean():.2f} max={obs_ctl.max():.0f} "
          f"| 误报显著 {n_sig_ctrl}/20（判 ≤1）→ {'OK' if ok2 else 'FAIL'}")
    out["g4b2"] = {"ctrl_mean_pairs": float(obs_ctl.mean()),
                   "ctrl_max_pairs": float(obs_ctl.max()),
                   "n_false_significant": int(n_sig_ctrl), "ok": ok2,
                   "note": "逐行 RA 置换（v1.2 修正: 块平移对库内配对无效）"}

    verdict = bool(ok1 and ok2)
    out["verdict"] = "PASS" if verdict else "FAIL"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g4b_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                         encoding="utf-8")
    print(f"G4b: {out['verdict']}  →  data/processed/upgrade/g4b_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())