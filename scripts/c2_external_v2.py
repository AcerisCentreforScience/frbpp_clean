"""Probe each candidate against the Blinkverse catalogue.

Usage: python scripts/c2_external_v2.py
Output: data/processed/upgrade/c2_candidates_v2.{parquet,csv}, c2_external_v2.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

BV = ROOT / "data" / "raw" / "blinkverse" / "blinkverse_sources.json"
V1 = ROOT / "data" / "processed" / "upgrade" / "c2_candidates_v1.parquet"
OUT = ROOT / "data" / "processed" / "upgrade"
SANITY = ["FRB20180916B", "FRB20121102A", "FRB20201124A"]


def haversine_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, [ra1, dec1, ra2, dec2])
    x = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    return float(np.degrees(np.arccos(x)) * 3600.0)


def main() -> int:
    bv = json.loads(BV.read_text(encoding="utf-8"))["data"]
    by_name = {}
    pos = []
    for x in bv:
        c = x.get("content", {})
        nm = str(c.get("source", "")).upper()
        rec = {"name": nm, "type": c.get("type"), "tel": c.get("telescope"),
               "ref": c.get("reference"), "dm": c.get("dm"),
               "ra": c.get("ra_degree"), "dec": c.get("dec_degree")}
        by_name[nm] = rec
        if np.isfinite(rec.get("ra") or np.nan) and np.isfinite(rec.get("dec") or np.nan):
            pos.append(rec)

    sanity = {s: by_name.get(s, {}).get("type") for s in SANITY}
    sanity_ok = all(v == "REPEAT_FIRST" for v in sanity.values())
    print(f"[锚点] 已知重复源: {sanity} → {'OK' if sanity_ok else 'FAIL'}")

    df = pd.read_parquet(V1)

    def probe(name, ra, dec):
        rec = by_name.get(str(name).upper())
        if rec is not None:
            return {"status": rec["type"], "tel": rec["tel"], "ref": rec["ref"],
                    "match": "name"}
        best, bd = None, np.inf
        if np.isfinite(ra) and np.isfinite(dec):
            for r in pos:
                d = haversine_arcsec(ra, dec, r["ra"], r["dec"])
                if d < bd:
                    best, bd = r, d
        if best is not None and bd <= 120.0:
            return {"status": f"近邻:{best['name']}({best['type']})",
                    "tel": best["tel"], "ref": best["ref"], "match": f"cone_{bd:.0f}\""}
        return {"status": "未收录", "tel": None, "ref": None, "match": None}

    rows = []
    for r in df.itertuples():
        p1 = probe(r.name1, r.ra1, r.dec1)
        p2 = probe(r.name2, r.ra2, r.dec2)
        flags = []
        for side, p, nm in (("1", p1, r.name1), ("2", p2, r.name2)):
            if p["status"] == "REPEAT_FIRST":
                flags.append(f"{nm}已为公开重复源")
            elif isinstance(p["status"], str) and p["status"].startswith("近邻:"):
                if "REPEAT_FIRST" in p["status"]:
                    flags.append(f"{nm}坐标近邻为公开重复源")
        if flags:
            verdict = "public_repeater:" + ";".join(flags)
        elif p1["status"] == "NON_REPEAT" and p2["status"] == "NON_REPEAT":
            verdict = "consistent_oneoff（双方均仍为非重复）"
        elif "未收录" in (p1["status"], p2["status"]):
            verdict = f"not_indexed（{r.name1}:{p1['status']}, {r.name2}:{p2['status']}）"
        else:
            verdict = f"mixed（{p1['status']}, {p2['status']}）"
        rows.append({
            "name1": r.name1, "name2": r.name2, "cls": r.cls, "dt_days": r.dt_days,
            "ddm": r.ddm, "ddm_sigma": r.ddm_sigma,
            "bv_status1": p1["status"], "bv_tel1": p1["tel"], "bv_ref1": p1["ref"],
            "bv_status2": p2["status"], "bv_tel2": p2["tel"], "bv_ref2": p2["ref"],
            "verdict": verdict,
        })
    v2 = pd.DataFrame(rows)
    v2.to_parquet(OUT / "c2_candidates_v2.parquet", index=False)
    v2.to_csv(OUT / "c2_candidates_v2.csv", index=False)

    n_pub = int(v2["verdict"].str.startswith("public_repeater").sum())
    n_ok = int(v2["verdict"].str.startswith("consistent_oneoff").sum())
    n_ni = int(v2["verdict"].str.startswith("not_indexed").sum())
    n_mx = len(v2) - n_pub - n_ok - n_ni
    print(f"[v2 汇总] 公开重复 {n_pub} | 一致(双非重复) {n_ok} | 未收录 {n_ni} | mixed {n_mx}")

    print("跨日强候选（逐对外部状态）:")
    for r in v2[v2["cls"] == "跨日强候选"].itertuples():
        print(f"  {r.name1} × {r.name2}: {r.verdict[:88]}")

    md = ["# C2 外部核查 v2（Blinkverse 在库比对，2026-10-08）", "",
          f"- Blinkverse 全源目录：4020 源（REPEAT_FIRST 101 / NON_REPEAT 3919），"
          f"来源 `data/raw/blinkverse/blinkverse_sources.json`（sha256 `3302A12F…`）",
          f"- 锚点：{SANITY} → REPEAT_FIRST {'OK' if sanity_ok else 'FAIL'}", "",
          f"- 28 对汇总：**公开重复 {n_pub}** | 一致双非重复 {n_ok} | 未收录 {n_ni} | mixed {n_mx}", ""]
    if n_pub or n_mx:
        md += ["## 需注意的对（公开重复 / mixed）", "",
               "| 对 | 类别 | bv1 | bv2 | 判定 |", "|---|---|---|---|---|"]
        for r in v2[v2["verdict"].str.startswith(("public_repeater", "mixed"))].itertuples():
            md.append(f"| {r.name1} × {r.name2} | {r.cls} | {r.bv_status1} | "
                      f"{r.bv_status2} | {r.verdict} |")
        md.append("")
    md += ["## 跨日强候选一览", "",
           "| 对 | Δt (d) | ΔDM/σ | bv1 | bv2 | 判定 |", "|---|---|---|---|---|---|"]
    for r in v2[v2["cls"] == "跨日强候选"].itertuples():
        s = f"{r.ddm_sigma:.1f}σ" if pd.notna(r.ddm_sigma) else "—"
        md.append(f"| {r.name1} × {r.name2} | {r.dt_days:.1f} | {s} | {r.bv_status1} | "
                  f"{r.bv_status2} | {r.verdict} |")
    (OUT / "c2_external_v2.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("→ c2_candidates_v2.parquet/csv + c2_external_v2.md")
    return 0 if sanity_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())