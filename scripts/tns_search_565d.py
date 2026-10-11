"""TNS Search API check for any-instrument detections inside the Cat2 blind window.

Usage: python scripts/tns_search_565d.py [--draws N] [--seed S] [--retries N] [--radius-arcsec R] [--out FILE]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config" / "tns.local.json"
OUTD = ROOT / "data" / "processed" / "upgrade"
LEDGER = OUTD / "tns_search_ledger.jsonl"
SUMMARY = OUTD / "tns_search_summary.md"
HOST = "www.wis-tns.org"
API = f"https://{HOST}/api/get/search"
API_OBJ = f"https://{HOST}/api/get/object"

PAIR = {
    "FRB20200511C": (20.187794, 51.151647),
    "FRB20211128E": (20.182712, 51.145837),
}
MID = (0.5 * (PAIR["FRB20200511C"][0] + PAIR["FRB20211128E"][0]),
       0.5 * (PAIR["FRB20200511C"][1] + PAIR["FRB20211128E"][1]))
DM_MID = 1834.95
BLIND = ("2023-09-15", "2026-09-05")
EPOCHS = {
    "k3_2025-01-02": ("2024-12-18", "2025-01-17"),
    "k4_2026-07-21": ("2026-07-06", "2026-08-05"),
}
RADII_ARCSEC = (600, 1800, 3600)

SURVIVING = [
    ("FRB20211127E", "FRB20220623D"), ("FRB20200919A", "FRB20210111A"),
    ("FRB20200814C", "FRB20210829A"), ("FRB20211116C", "FRB20211129A"),
    ("FRB20220607C", "FRB20220621D"), ("FRB20200220I", "FRB20230721C"),
]

CTRL_NAMES_POS = ["FRB20180916B", "FRB20220912A", "FRB20201124A"]
CTRL_NAMES_NEG = ["FRB20990101Z", "FRB20881231Q"]
BLANK_FIELD = (198.0, -60.0)


def marker(bot_id, bot_name) -> str:
    return f'tns_marker{{"tns_id": "{bot_id}", "type": "bot", "name": "{bot_name}"}}'


def sexa(ra_deg, dec_deg):
    c = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    return c.ra.to_string(unit=u.hour, sep=":", pad=True, precision=2), \
        c.dec.to_string(unit=u.deg, sep=":", pad=True, precision=1, alwayssign=True)


class TNS:
    def __init__(self, bot_id, bot_name, api_key, sleep=1.0):
        self.h = {"User-Agent": marker(bot_id, bot_name)}
        self.k = api_key
        self.sleep = sleep
        self.n_calls = 0

    def search(self, **data):
        payload = {"api_key": self.k, "data": json.dumps(data)}
        r = requests.post(API, headers=self.h, data=payload, timeout=60)
        self.n_calls += 1
        try:
            j = r.json()
        except Exception:
            return r.status_code, None, {"_raw": r.text[:300]}, {}
        code = j.get("id_code")
        rl = {k: v for k, v in r.headers.items()
              if k.lower().startswith(("x-rate-limit", "x-cone-rate-limit"))}
        if code == 401:
            return code, None, j, rl
        if code == 403:
            reset = self._reset(r.headers)
            if reset:
                time.sleep(min(reset + 1, 120))
            return code, None, j, rl
        d = j.get("data")
        objs = d.get("reply") or [] if isinstance(d, dict) else (d or [])
        names = []
        for o in objs:
            if isinstance(o, dict):
                names.append(f"{o.get('prefix', '')}{o.get('objname', '')}")
            elif o:
                names.append(str(o))
        rem = [v for k, v in r.headers.items() if k.lower().endswith("rate-limit-remaining")]
        if any(x.isdigit() and int(x) <= 2 for x in rem):
            reset = self._reset(r.headers)
            if reset:
                time.sleep(min(reset + 1, 120))
        time.sleep(self.sleep)
        return code, names, j, rl

    @staticmethod
    def _reset(h):
        for k in ("x-rate-limit-reset", "x-cone-rate-limit-reset"):
            v = h.get(k)
            if v and str(v).isdigit() and int(v) > 0:
                return int(v)
        return None


    def get_object(self, objname):
        payload = {"api_key": self.k,
                   "data": json.dumps({"objname": objname, "photometry": "0", "spectra": "0"})}
        r = requests.post(API_OBJ, headers=self.h, data=payload, timeout=60)
        self.n_calls += 1
        try:
            j = r.json()
        except Exception:
            return r.status_code, None
        code = j.get("id_code")
        d = j.get("data")
        if isinstance(d, list):
            d = d[0] if d else None
        obj = d.get("reply", d) if isinstance(d, dict) else d
        rem = [v for k, v in r.headers.items() if k.lower().endswith("rate-limit-remaining")]
        if any(x.isdigit() and int(x) <= 2 for x in rem):
            reset = self._reset(r.headers)
            if reset:
                time.sleep(min(reset + 1, 120))
        time.sleep(self.sleep)
        return code, obj


def build_queries(ra=None, dec=None, radius=None, units="arcsec",
                  dates=None, dm=None, **kw):
    q = {"include_frb": "1", "public_search": "1"}
    if ra is not None:
        sra, sdec = sexa(ra, dec)
        q.update({"ra": sra, "dec": sdec, "radius": str(radius), "units": units})
    if dates:
        q.update({"discovery_date_start": dates[0], "discovery_date_end": dates[1]})
    if dm:
        q.update({"frb_dm_range_min": str(dm[0]), "frb_dm_range_max": str(dm[1])})
    q.update(kw)
    return q


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="只打印将发出的查询，不联网")
    a = ap.parse_args()
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    tns = TNS(cfg.get("tns_id", 0), cfg["bot_name"], cfg["api_key"])

    plan = []
    for n in CTRL_NAMES_POS:
        plan.append(("CTRL_pos_name", build_queries(objname=n, objname_exact_match=0)))
    for n in CTRL_NAMES_NEG:
        plan.append(("CTRL_neg_name", build_queries(objname=n)))
    plan.append(("CTRL_neg_blank", build_queries(
        ra=BLANK_FIELD[0], dec=BLANK_FIELD[1], radius=3600)))
    plan.append(("CTRL_pos_cone_180916B", build_queries(
        ra=29.503, dec=65.717, radius=3600)))

    for nm, (ra, dec) in PAIR.items():
        for R in RADII_ARCSEC:
            plan.append((f"cone_{nm}_{R}", build_queries(ra=ra, dec=dec, radius=R)))
    for R in RADII_ARCSEC:
        plan.append((f"cone_mid_{R}_blind", build_queries(
            ra=MID[0], dec=MID[1], radius=R, dates=BLIND)))
    for lab, win in EPOCHS.items():
        plan.append((f"epoch_{lab}", build_queries(
            ra=MID[0], dec=MID[1], radius=3600, dates=win)))
    plan.append(("dm_fingerprint", build_queries(
        frb_dm_range_min=DM_MID - 20, frb_dm_range_max=DM_MID + 20)))
    for nm in PAIR:
        plan.append((f"tns_entry_{nm}", build_queries(objname=nm, objname_exact_match=1)))

    for a1, a2 in SURVIVING:
        for nm in (a1, a2):
            plan.append((f"surv_{nm}", build_queries(objname=nm)))

    if a.dry_run:
        print(f"[dry-run] 计划查询 {len(plan)} 条；端点 {API}")
        for lab, q in plan[:12]:
            print(f"  {lab:26s} {json.dumps(q, ensure_ascii=False)}")
        print("  ...")
        return 0

    OUTD.mkdir(parents=True, exist_ok=True)
    rows, auth_failed = [], False
    with LEDGER.open("w", encoding="utf-8") as fh:
        for lab, q in plan:
            code, names, raw, rl = tns.search(**q)
            if code == 401:
                print(f"[{lab}] id_code=401 Unauthorized —— **凭证被拒，停止检索**")
                print("        未回显任何凭证。请核对 bot 的 API Key（MY ACCOUNT → API keys，"
                      "须是 FRB11 这个 bot 的 key，不是用户 key；生产与沙盒账号相互独立）")
                auth_failed = True
                break
            rec = {"ts": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                   "label": lab, "query": q, "id_code": code,
                   "n": (len(names) if names is not None else None),
                   "names": (names or [])[:60], "rate": rl}
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fh.flush()
            rows.append(rec)
            shown = ",".join((rec["names"] or [])[:8])
            print(f"[{code}] {lab:28s} n={str(rec['n']):>4s}  {shown}")
            if code == 400:
                print(f"      400 反馈: {json.dumps(raw, ensure_ascii=False)[:220]}")

    if auth_failed:
        return 2

    pos = [r for r in rows if r["label"].startswith("CTRL_pos")]
    neg = [r for r in rows if r["label"].startswith("CTRL_neg")]
    power = any((r["n"] or 0) > 0 for r in pos)
    clean_neg = all((r["n"] or 0) == 0 for r in neg)
    main_hits = [r for r in rows if not r["label"].startswith("CTRL")]
    found = [r for r in main_hits if (r["n"] or 0) > 0]

    lines = [f"# TNS 检索汇总（{dt.date.today()}）", "",
             f"- 端点 `{API}`；调用 {tns.n_calls} 次；查询 {len(rows)} 条",
             f"- **判据功效**：正对照 {'有命中 ✓' if power else '全空 ✗ → 本轮所有 0 命中不得作为证据（未验证）'}",
             f"- 负对照干净：{'是 ✓' if clean_neg else '否 ✗ → 圆锥/参数过松，需收紧'}",
             f"- 主查命中条目：{len(found)}/{len(main_hits)}", ""]
    for r in main_hits:
        if r["n"]:
            lines.append(f"  - **{r['label']}** n={r['n']}: {','.join(r['names'][:12])}")
    if not found:
        lines.append("  （主查全部为空 —— 仅当上面功效 ✓ 时才是有意义的负结果）")
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if power else 3


if __name__ == "__main__":
    raise SystemExit(main())
