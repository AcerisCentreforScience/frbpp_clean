"""Audit of the event-table de-duplication conventions quoted in Section 2.

Usage: python scripts/event_table_dedup_audit.py
Output: data/processed/upgrade/event_table_audit.json
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
OUT = ROOT / "data" / "processed" / "upgrade"

SAME_EVENT_DT_DAYS = 0.05


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _clean_rep(s) -> np.ndarray:
    x = s.astype("string").str.strip().fillna("")
    return x.mask(x.isin(["-9999", "nan", "NaN", "None", "none", ""]), "").to_numpy(dtype=object)


def main() -> dict:
    c1 = pd.read_parquet(INTERIM / "cat1_bursts.parquet")
    c2 = pd.read_parquet(INTERIM / "cat2_bursts.parquet")

    m1, m2 = _num(c1.mjd_400).to_numpy(), _num(c2.mjd_400).to_numpy()
    ok1 = (_num(c1.ra).notna() & _num(c1.dec).notna() & np.isfinite(m1)).to_numpy()
    ok2 = (_num(c2.ra).notna() & _num(c2.dec).notna() & np.isfinite(m2)).to_numpy()

    c2o = c2[ok2].reset_index(drop=True)
    m2o = _num(c2o.mjd_400).to_numpy()
    dm2o = _num(c2o.dm_fitb).to_numpy()
    n2 = c2o.tns_name.astype(str).to_numpy()
    n1 = c1.tns_name.astype(str).to_numpy()
    dm1 = _num(c1.dm_fitb).to_numpy()
    sub1 = _num(c1.sub_num).fillna(0).to_numpy()

    idx: dict[str, list[int]] = {}
    for k, nm in enumerate(n2):
        idx.setdefault(nm, []).append(k)

    dup1 = np.zeros(len(c1), dtype=bool)
    for r in range(len(c1)):
        if not ok1[r]:
            continue
        for k in idx.get(n1[r], []):
            if abs(m1[r] - m2o[k]) < SAME_EVENT_DT_DAYS:
                dup1[r] = True
                break

    keep1 = ok1 & ~dup1
    r1, r2 = _clean_rep(c1.repeater_name), _clean_rep(c2o.repeater_name)

    row_all = int(ok1.sum()) + int(ok2.sum())
    row_rep = int((r1 != "").sum()) + int((r2 != "").sum())
    row_off = row_all - row_rep
    oo_row = np.concatenate([r1 == "", r2 == ""])
    dm_ok_row = np.concatenate([np.isfinite(dm1), np.isfinite(dm2o)])
    row_off_finite = int(np.sum(oo_row & dm_ok_row))

    all_d = int(ok2.sum()) + int(keep1.sum())
    rep_d = int((r2 != "").sum()) + int(np.sum(keep1 & (r1 != "")))
    off_d = all_d - rep_d
    oo_d = np.concatenate([(r1 == "") & keep1, (r2 == "")])
    dm_ok_d = np.concatenate([np.isfinite(dm1) & keep1, np.isfinite(dm2o)])
    off_d_finite = int(np.sum(oo_d & dm_ok_d))

    out = {
        "built_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "criteria": "同名 且 |Δt|<0.05 d = 同一暴（与 crossmatch.SAME_EVENT_DT_DAYS 一致）；保留 Cat2 为权威行",
        "note": (
            "build_event_table 只在目录行级 concat，去重铁律在配对层生效。"
            "故 row-level 计数含 Cat1/Cat2 重复收录行；dedup 计数为独立物理事件数。"
        ),
        "row_level": {
            "cat1_rows": int(len(c1)),
            "cat2_rows_ok": int(ok2.sum()),
            "cat2_rows_dropped_no_coords": int((~ok2).sum()),
            "all": row_all,
            "repeater": row_rep,
            "one_off": row_off,
            "one_off_finite_dm": row_off_finite,
        },
        "dedup": {
            "cat1_rows_matched_in_cat2": int(dup1.sum()),
            "cat1_rows_only_in_cat1": int(keep1.sum()),
            "cat1_sub_burst_rows_total": int((sub1 >= 1).sum()),
            "cat1_sub_burst_rows_matched": int(np.sum((sub1 >= 1) & dup1)),
            "all": all_d,
            "repeater": rep_d,
            "one_off": off_d,
            "one_off_finite_dm": off_d_finite,
        },
        "inflation": {
            "all_rows": row_all - all_d,
            "one_off_rows": row_off - off_d,
            "all_frac": (row_all - all_d) / all_d,
        },
        "impact": (
            "候选对/假设计数不受影响——重复行在配对层已被 same_event_mask 剔除；"
            "受影响的是所有'事件池计数'陈述与 §4.7 误差卷积 MC 的抽样池。"
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / "event_table_audit.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print(f"\n-> {dst}")
    return out


if __name__ == "__main__":
    main()
