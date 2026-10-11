"""Literature check for published RM or scattering measurements of the pair.

Usage: python scripts/rm_ads_565d.py
Output: data/processed/upgrade/pair565_rm_ads.json
"""
from __future__ import annotations

import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "config" / "ads.local.json"
OUT = ROOT / "data" / "processed" / "upgrade"
URL = "https://api.adsabs.harvard.edu/v1/search/query"

QUERIES = [
    ("20200511C_RM", 'full:"FRB20200511C" AND full:"rotation measure"'),
    ("20211128E_RM", 'full:"FRB20211128E" AND full:"rotation measure"'),
    ("20200511C_scat", 'full:"FRB20200511C" AND full:"scattering"'),
    ("20211128E_scat", 'full:"FRB20211128E" AND full:"scattering"'),
    ("20200511C_pol", 'full:"FRB20200511C" AND full:"polarization"'),
    ("CTRL_20201124A_RM", 'full:"FRB20201124A" AND full:"rotation measure"'),
    ("CTRL_180916_RM", 'full:"FRB20180916B" AND full:"rotation measure"'),
]


def main() -> int:
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    hdr = {"Authorization": f"Bearer {cfg['api_token']}"}
    out = {}
    for tag, q in QUERIES:
        r = requests.get(URL, headers=hdr, params={"q": q, "fl": "title,year,bibcode", "rows": 5}, timeout=40)
        j = r.json()["response"]
        out[tag] = {"q": q, "numFound": j["numFound"], "docs": j["docs"][:3]}
        print(f"{tag:22s} numFound={j['numFound']:3d} {[d.get('bibcode') for d in j['docs'][:3]]}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair565_rm_ads.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("-> data/processed/upgrade/pair565_rm_ads.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
