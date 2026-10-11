"""Download third-party FRB catalogues via VizieR TAP into data/raw.

Usage: python scripts/fetch_external_catalogs.py
"""
import hashlib
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"


def fetch(table: str) -> dict:
    r = requests.get(TAP, params={"request": "doQuery", "lang": "ADQL",
                                  "format": "json", "query": f'SELECT * FROM "{table}"'},
                     timeout=180)
    r.raise_for_status()
    j = r.json()
    cols = [m["name"] for m in j["metadata"]]
    rows = [dict(zip(cols, row)) for row in j["data"]]
    return {"table": table, "columns": cols, "n_rows": len(rows), "rows": rows}


JOBS = [
    ("J/A+A/693/A279/frbs", ROOT / "data" / "raw" / "apertif", "apertif_frbs.json"),
    ("J/ApJ/955/142/table1", ROOT / "data" / "raw" / "fast", "fast_20220912A_bursts.json"),
]

for table, folder, fname in JOBS:
    d = fetch(table)
    folder.mkdir(parents=True, exist_ok=True)
    body = json.dumps(d, ensure_ascii=False, indent=1)
    p = folder / fname
    p.write_text(body, encoding="utf-8")
    sha = hashlib.sha256(body.encode()).hexdigest()
    print(f"[{table}] rows={d['n_rows']} → {p.relative_to(ROOT)}  sha256={sha[:16]}…")