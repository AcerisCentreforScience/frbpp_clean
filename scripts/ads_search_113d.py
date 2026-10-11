"""Literature check for the 113.69 d pair.

Usage: python scripts/ads_search_113d.py
"""
from __future__ import annotations

import datetime as dt
import json
import re
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config" / "ads.local.json"
OUT = ROOT / "data" / "processed" / "upgrade"
URL = "https://api.adsabs.harvard.edu/v1/search/query"

QUERIES = [
    ("name_20200919A", '"FRB20200919A"'),
    ("name_20210111A", '"FRB20210111A"'),
    ("name_spaced", '"FRB 20200919A"'),
    ("short_id", '"20200919A"'),
    ("fulltext_1", 'full:"FRB20200919A"'),
    ("fulltext_2", 'full:"FRB20210111A"'),
    ("coords_pair", 'full:"FRB20200919A" OR full:"FRB20210111A"'),
    ("CTRL_known_repeater", '"FRB20180916B"'),
    ("CTRL_known_repeater2", '"FRB20220912A"'),
    ("CTRL_nonsense", '"FRB20990101Z"'),
    ("CTRL_nonsense2", '"FRB20991231Q"'),
    ("CTRL_bibcode_field", 'bibstem:"GCN" AND full:"FRB20200919A"'),
]


def main() -> int:
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    tok = cfg["api_token"]
    hdr = {"Authorization": f"Bearer {tok}"}
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "queries": {}}
    print(f"查询 {len(QUERIES)} 组 …")
    for tag, q in QUERIES:
        try:
            r = requests.get(URL, headers=hdr,
                             params={"q": q, "fl": "title,year,bibcode,first_author",
                                     "rows": 25}, timeout=40)
            if r.status_code != 200:
                out["queries"][tag] = {"q": q, "status": r.status_code,
                                       "body": r.text[:200]}
                print(f"  [{tag}] {r.status_code} {r.text[:120]}")
                continue
            j = r.json()
            docs = j["response"]["docs"]
            out["queries"][tag] = {"q": q, "numFound": j["response"]["numFound"],
                                   "docs": docs}
            print(f"  [{tag}] numFound={j['response']['numFound']}")
            for d in docs[:6]:
                print(f"       {d.get('year')} | {d.get('first_author')} | "
                      f"{str(d.get('title'))[:100]}")
        except Exception as e:  # noqa: BLE001
            out["queries"][tag] = {"q": q, "error": str(e)[:200]}
            print(f"  [{tag}] ERROR {str(e)[:120]}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair113_ads.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\n→ data/processed/upgrade/pair113_ads.json")

    print("\n===== 命中详情 =====")
    for probe_tag in ("name_20200919A", "name_20210111A"):
        hits = out["queries"].get(probe_tag, {}).get("docs", [])
        for d in hits:
            bib = d.get("bibcode")
            r = requests.get(URL, headers=hdr,
                             params={"q": f'bibcode:"{bib}"',
                                     "fl": "title,abstract,year,author,bibcode,body,esources,keyword",
                                     "rows": 1}, timeout=40)
            if r.status_code != 200:
                print(f"  [{bib}] detail {r.status_code}")
                continue
            dd = r.json()["response"]["docs"][0]
            body = str(dd.get("body", ""))
            print(f"  bibcode={bib}\n  title={dd.get('title')}\n  year={dd.get('year')}")
            print(f"  abstract={str(dd.get('abstract'))[:400]}")
            print(f"  body_len={len(body)} esources={dd.get('esources')}")
            for term in ("20200919A", "20210111A"):
                i = body.find(term)
                print(f"    body 含 {term}? pos={i}" +
                      (f"  …{body[max(0,i-260):i+120]}…" if i >= 0 else ""))
            kw = dd.get("keyword") or []
            frb_names = [k for k in kw if re.fullmatch(r"FRB\d{8}[A-Z]", str(k))]
            yrs = sorted({str(n)[3:7] for n in frb_names})
            print(f"  keyword 总数={len(kw)}  其中 FRB 名={len(frb_names)}  年份跨度={yrs}")
            if len(frb_names) > 50:
                print("  ⚠ 判定：**批量登记工件** —— 该命中由 TNS 批量条目携带的名字清单触发，"
                      "不代表有再探测")

    print("\n===== 特异性检验（同一 TNS 报告是否吞下所有 CHIME 源名）=====")
    probes = ["FRB20190303A", "FRB20201124A", "FRB20190117A", "FRB20180814A",
              "FRB20190417A", "FRB20200120E", "FRB20190915D", "FRB20190212A",
              "FRB20180916B", "FRB20220912A",
              "FRB20200919A", "FRB20210111A", "FRB20200511C", "FRB20211128E",
              "FRB20991231Q", "FRB20880101Z"]
    for nm in probes:
        rr = requests.get(URL, headers=hdr,
                          params={"q": f'"{nm}"', "fl": "bibcode,title,year",
                                  "rows": 3}, timeout=30)
        if rr.status_code != 200:
            print(f"  {nm:15s} -> HTTP {rr.status_code}")
            continue
        jj = rr.json()["response"]
        bibs = [d.get("bibcode") for d in jj["docs"]][:3]
        print(f"  {nm:15s} numFound={jj['numFound']:4d}  {bibs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
