"""Verify each bibitem against ADS records.
Superseded by verify_refs_ads2.py.
"""
import json
import os
import re
import sys
import time
import pathlib

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
os.environ.setdefault("NO_PROXY", "*")
TEX = ROOT / "docs/submission/raa_manuscript.tex"
TOKEN = json.loads((ROOT / "config/ads.local.json").read_text(encoding="utf-8"))["api_token"]
H = {"Authorization": f"Bearer {TOKEN}", "User-Agent": "frbpp-refcheck/1.0"}
URL = "https://api.adsabs.harvard.edu/v1/search/query"
FL = "bibcode,year,title,v,pub,pubdate,doctypes,esourceid"

lines = TEX.read_text(encoding="utf-8").splitlines()
entries = []
i = 0
while i < len(lines):
    m = re.match(r"\s*\\bibitem\{([^}]+)\}", lines[i])
    if m:
        key, buf, j = m.group(1), [re.sub(r"^\s*\\bibitem\{[^}]+\}\s*", "", lines[i])], i + 1
        while j < len(lines) and not re.match(r"\s*\\bibitem\{", lines[j]) and not re.match(
                r"\s*\\end\{thebibliography\}", lines[j]):
            buf.append(lines[j].strip())
            j += 1
            if len(buf) > 7:
                break
        entries.append((key, " ".join(b for b in buf if b)))
        i = j
    else:
        i += 1

print(f"解析到 {len(entries)} 条\n")


def parse(txt):
    t = txt.replace("\\&", "&")
    t = re.sub(r"\\(cite[tp]?\{[^}]*\}|textbf|emph|url|doi)\{?[^}]*\}?", " ", t)
    t = re.sub(r"[{}]", "", t).replace("  ", " ").strip()
    year = re.search(r"\b(19[5-9]\d|20[0-4]\d)\b", t)
    yr = year.group(1) if year else None
    vp = re.search(r",\s*(19|20)\d\d,\s*([A-Za-z\-\. ]{1,28}?)\s*(\d+)\s*:\s*([0-9LA-Za-z\-\u2013]+)", t)
    journal = None
    vol = page = None
    if vp:
        journal, vol, page = vp.group(2).strip(), vp.group(3), vp.group(4).split("-")[0].split("\u2013")[0]
    else:
        jm = re.search(r",\s*(?:19|20)\d\d\s*$", t)
    first = t.split(",")[0].strip()
    arx = re.search(r"arXiv:([0-9a-zA-Z\.\-/]+)", t)
    return dict(year=yr, journal=journal, vol=vol, page=page, first_author=first,
                arxiv=arx.group(1) if arx else None, raw=t)


def ads_q(q):
    r = requests.get(URL, headers=H, params={"q": q, "fl": FL, "rows": 3}, timeout=45)
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    return r.json().get("response", {}).get("documents", []), None


problems, ok = [], 0
for key, txt in entries:
    p = parse(txt)
    doc = None
    if p["arxiv"]:
        docs, err = ads_q(f'property:{"arxiv"} AND abs:"{p["arxiv"]}"')
        if docs is None:
            print(f"[{key}] ADS 查询失败 {err}")
            problems.append((key, "ADS-QUERY-FAIL", err))
            continue
        if docs:
            doc = docs[0]
    if doc is None:
        title_bits = re.findall(r"[A-Z][^.]{28,}", p["raw"])[:1]
        term = title_bits[0][:70].strip() if title_bits else p["raw"][:60]
        docs, err = ads_q(f'title:"{term}"')
        if docs is None:
            print(f"[{key}] ADS 查询失败 {err}")
            problems.append((key, "ADS-QUERY-FAIL", err))
            continue
        if docs:
            doc = max(docs, key=lambda d: 1 if d.get("esourceid") else 0)
    if doc is None:
        print(f"[{key}]  NOT-FOUND  稿内: {p['raw'][:110]}")
        problems.append((key, "NOT-FOUND-IN-ADS", p["raw"][:90]))
        continue

    ay = str(doc.get("year") or "")
    aj = (doc.get("pub") or "")
    av = str(doc.get("v") or "")
    ap = str(doc.get("page") or "")
    notes = []
    if p["year"] and ay and p["year"] != ay:
        notes.append(f"年份 稿{p['year']} vs ADS {ay}")
    if p["vol"] and av and p["vol"] != av:
        notes.append(f"卷号 稿{p['vol']} vs ADS {av}")
    if p["page"] and ap and not ap.startswith(p["page"]) and p["page"] not in ap:
        notes.append(f"页码 稿{p['page']} vs ADS {ap}")
    tag = "OK " if not notes else "CHK"
    if notes:
        problems.append((key, "FIELD-MISMATCH", "; ".join(notes)))
    else:
        ok += 1
    print(f"[{key}] {tag} ADS={doc.get('bibcode')} {aj[:34]:34s} v{av} p{ap} | "
          f"{' ; '.join(notes) if notes else '一致'}")
    time.sleep(0.4)

print(f"\n完全一致 {ok}/{len(entries)}")
print("需人工定夺:")
for k, kind, detail in problems:
    print(f"  - {k:18s} {kind:18s} {str(detail)[:120]}")
